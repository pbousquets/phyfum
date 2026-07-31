# Crypt-division model averaging

Created: 2026-07-28

Last updated: 2026-07-31

## Purpose and model definition

PHYFUM can sample the crypt-division model jointly with its other MCMC
parameters instead of running four independent analyses and comparing them
afterward. One global scalar `divisionModel` parameter controls the pruning
calculation at every internal node:

| Value | Model | Pruning implementation |
|---:|---|---|
| 0 | identity | `GeneralCenancestorLikelihoodCore` |
| 1 | budding | `BuddingCenancestorLikelihoodCore` |
| 2 | fission | `FissionCenancestorLikelihoodCore` |
| 3 | split fission | `SplitFissionCenancestorLikelihoodCore` |

The selector is global, not node-specific. The implementation is Java-only and
currently accepts `AFsequence` data. Existing fixed-model
`cenancestorTreeLikelihood` XML remains compatible and unchanged.

## Source map

| Responsibility | File |
|---|---|
| Tree-likelihood model and selector events | `src/dr/evomodel/treelikelihood/ModelAveragingCenancestorTreeLikelihood.java` |
| Shared-buffer likelihood core and pruning dispatch | `src/dr/evomodel/treelikelihood/ModelAveragingCenancestorLikelihoodCore.java` |
| Categorical prior | `src/dr/evomodel/treelikelihood/DivisionModelPrior.java` |
| Tree-likelihood XML parser | `src/dr/evomodelxml/treelikelihood/ModelAveragingCenancestorTreeLikelihoodParser.java` |
| Prior XML parser | `src/dr/evomodelxml/treelikelihood/DivisionModelPriorParser.java` |
| Parser registration | `src/dr/app/beast/release_parsers.properties` |
| Runnable example | `examples/release/flipflop/modelAveraging.xml` |
| Core equivalence and state tests | `src/test/dr/evomodel/treelikelihood/ModelAveragingCenancestorLikelihoodCoreTest.java` |
| Complete likelihood tests | `src/test/dr/evomodel/treelikelihood/ModelAveragingCenancestorTreeLikelihoodTest.java` |
| Prior and operator tests | `src/test/dr/evomodel/treelikelihood/DivisionModelPriorTest.java` |

## Architecture and buffer ownership

`ModelAveragingCenancestorTreeLikelihood` extends `CenancestorTreeLikelihood`
and supplies a `ModelAveragingCenancestorLikelihoodCore`. The model-averaging
core extends `GeneralCenancestorLikelihoodCore`, so the existing traversal,
cenancestor handling, scaling, and double-buffer state machinery remain in use.

There is exactly one set of tree-sized likelihood-core storage:

- matrices;
- partials;
- tip states;
- scaling factors;
- current and stored matrix/partial buffer indices.

These arrays belong to `ModelAveragingCenancestorLikelihoodCore` and are
allocated by its call to `super.initialize(...)`.

The four objects in `divisionCores` are algorithm delegates only. Their
`initialize(...)` methods are deliberately not called, so their inherited
matrices, partials, states, scaling buffers, and current/stored buffer-index
arrays remain `null`. They receive the model-averaging core's arrays as
arguments to their protected pruning methods. Only their loop dimensions and
model-specific lookup/scratch data are initialized.

The delegates are likelihood-core objects, not BEAST `Model` objects. They are
not registered with `addModel(...)` or `addVariable(...)`, have no dirty flags,
and receive no model or variable listener events. A selector change marks the
single tree likelihood's node flags; it does not mark four independent cores.
During traversal, only the cached selected delegate executes.

Model-specific lookup tables are initialized once at setup so every model is
ready when first proposed. This includes budding probabilities and the fission
and split-fission combination tables. These tables are distinct biological
data and would still be required by a fused implementation; they are not
duplicate tree likelihood buffers.

The pruning path is:

```text
ModelAveragingCenancestorLikelihoodCore matrices/partials
        -> inherited calculatePartials(...)
        -> overridden pruning dispatch
        -> cached selected delegate's protected pruning method
```

Delegation preserves the fixed cores' numerical calculations and their current
`UnsupportedOperationException` behavior for unsupported pruning variants.

## Selector validation and cached dispatch

The constructor requires a scalar selector, validates that its current value
is an exact integer from 0 through 3, and attaches bounds `[0,3]`. Bounds alone
do not enforce integrality because BEAST's `Parameter` stores doubles; the
configured integer operator supplies valid MCMC proposals.

`getDivisionModelIndex(...)` both validates and returns the selector.
`validateDivisionModel(...)` is used where only validation is intended, making
discarded-return call sites explicit.

Pruning does not read the parameter for every node. The core caches:

```text
currentDivisionModel  selected core for the proposed/current state
storedDivisionModel   selected core at the last MCMC store point
```

`getDivisionCore()` indexes `divisionCores` with `currentDivisionModel`, so
node-level pruning performs only a cached array lookup.

## MCMC proposal, acceptance, and rejection lifecycle

The order of operations is important:

1. Before an operator runs, BEAST calls `storeModelState()`.
2. Core `storeState()` stores matrix/partial buffer indices and copies
   `currentDivisionModel` into `storedDivisionModel`.
3. `uniformIntegerOperator` calls `Parameter.setParameterValue(...)`.
4. The parameter writes the proposed value before firing its listener.
5. `isUpdatingDivisionModel()` therefore reads the new proposal with
   `getDivisionModelIndex(...)`, while `currentDivisionModel` still contains
   the pre-proposal model.
6. If the two models differ, the cache changes to the proposal and every tree
   node is marked for recalculation under that model.

If the proposal is accepted, `currentDivisionModel` remains the selected model.
The next MCMC store operation refreshes `storedDivisionModel`.

If the proposal is rejected, BEAST restores the parameter before restoring the
tree likelihood's additional state. The core then restores
`currentDivisionModel = storedDivisionModel` together with its stored
matrix/partial buffer indices. A forced post-restore pruning test verifies that
the restored cache selects the original model, rather than merely observing
previously stored partial values.

`CenancestorTreeLikelihood` normally invokes core `storeState()` and
`restoreState()` only when `storePartials=true`. When `storePartials=false`, the
model-averaging tree likelihood explicitly stores and restores the selector
cache; the parent likelihood marks nodes for recalculation as usual. Both modes
are covered by rejection tests.

## Avoiding unnecessary self-proposal recalculation

BEAST's existing `uniformIntegerOperator` samples uniformly over all four
categories, including the current value. This is a symmetric proposal with log
Hastings ratio `0.0`.

`isUpdatingDivisionModel()` compares the proposed selector with the cached
current selector. If they are equal, it returns `false` without changing the
cache. The tree-likelihood handler then leaves all node flags unchanged, so the
MCMC iteration represents the same state without repeating a full-tree pruning
calculation. The proposal distribution and Hastings ratio are unchanged.

## XML interface

The new likelihood is declared with:

```xml
<modelAveragingCenancestorTreeLikelihood id="treeLikelihood">
    ...
    <divisionModel>
        <parameter id="divisionModel" value="0" lower="0" upper="3"/>
    </divisionModel>
    ...
</modelAveragingCenancestorTreeLikelihood>
```

The selector is operated on with BEAST's existing uniform integer operator:

```xml
<uniformIntegerOperator lower="0" upper="3" weight="3.0">
    <parameter idref="divisionModel"/>
</uniformIntegerOperator>
```

The selector should be included in an output log so posterior model
probabilities can be estimated from its sampled frequencies.

## Division-model prior

The categorical prior is deliberately separate from the tree likelihood so
tree-likelihood values remain comparable with fixed-model runs. It must be
placed explicitly in the MCMC `<prior>` section:

```xml
<prior id="prior">
    ...
    <divisionModelPrior weights="1.0 1.0 1.0 1.0">
        <parameter idref="divisionModel"/>
    </divisionModelPrior>
</prior>
```

Weights are fixed relative weights in identity, budding, fission, and split
fission order and are normalized internally. Omitting `weights` gives equal
prior probabilities. A zero weight excludes a model. Negative, non-finite,
wrong-length, and all-zero vectors are rejected.

## Maintenance invariants

Future changes should preserve the following:

- Do not call `initialize(...)` on objects in `divisionCores`; doing so would
  allocate redundant tree-sized buffers.
- Keep `divisionCores` in exactly the same order as the public numeric mapping
  and prior weights.
- Update `currentDivisionModel` once in the selector listener, never during
  every node-level pruning call.
- Store and restore the cached selector for both `storePartials` modes.
- A real selector change must invalidate all nodes because shared partials were
  calculated under the previous model.
- A same-selector proposal must not invalidate nodes.
- Preserve fixed-core behavior, including unsupported pruning variants.
- Adding another model requires updating `MODEL_COUNT`, numeric constants,
  delegate construction, prior validation/documentation, XML comments, and
  equivalence tests.

## Verification

The following checks passed:

- `ant compile-all`.
- Four core tests covering all fixed-core equivalence comparisons, cached
  update detection, invalid selectors, and forced pruning after cache restore.
- Two complete tree-likelihood tests covering all four fixed likelihoods and
  proposal rejection with both `storePartials=true` and `false`.
- Five prior/operator tests covering equal and unequal weights, zero weights,
  invalid weights, integer support, and the zero Hastings ratio.
- An initial 10,000-state smoke run with:

  ```text
  dr.app.beast.BeastMain -beagle_off -overwrite -seed 12345 \
      examples/release/flipflop/modelAveraging.xml
  ```

  The run completed in 0.296 seconds. All four selector values appeared in the
  log. The division-model operator made 9,900 proposals with acceptance
  probability 0.7921.

The existing `ant junit_flipflop` target retained its baseline result: all
discoverable tests pass, while `TestSubstitutionModelEmpiricalFrequencies`
cannot be discovered because the source declares package
`dr.evomodel.flipflop` but the Ant target requests
`test.dr.evomodel.flipflop`. This pre-existing unrelated test was not changed.

## Deferred validation

Posterior model frequencies have not yet been compared with the existing
independent fixed-model simulation runs. That is the next validation phase.
