# Crypt-division model averaging

Date: 2026-07-28

## Implementation

PHYFUM can now sample one crypt-division model within the MCMC. The scalar
`divisionModel` parameter applies to every internal node and uses this mapping:

| Value | Model |
|---:|---|
| 0 | identity |
| 1 | budding |
| 2 | fission |
| 3 | split fission |

`ModelAveragingCenancestorTreeLikelihood` extends the existing cenancestor tree
likelihood. Its core owns one set of partial and matrix buffers and delegates
the pruning calculations to the selected existing likelihood core. A selector
change marks every node for recalculation. BEAST's normal parameter and
likelihood-core store/restore mechanisms therefore restore the previous
partials after a rejected model proposal.

The existing `cenancestorTreeLikelihood` element and its fixed
`divisionModel` attribute are unchanged.

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

This proposal is symmetric over the four categories, including a possible
self-transition, so its log Hastings ratio is `0.0`.

The categorical model prior is a separate component and belongs explicitly in
the MCMC `<prior>` section:

```xml
<prior id="prior">
    ...
    <divisionModelPrior weights="1.0 1.0 1.0 1.0">
        <parameter idref="divisionModel"/>
    </divisionModelPrior>
</prior>
```

The values are fixed relative weights in identity, budding, fission, and split
fission order. They are normalized internally. Omitting `weights` gives equal
prior probability. A zero weight excludes a model; negative, non-finite,
wrong-length, and all-zero weight vectors are rejected.

A compact runnable configuration is provided in
`examples/release/flipflop/modelAveraging.xml`.

## Verification

The following checks passed:

- Clean Java compilation with `ant compile-all`.
- Three core tests comparing every dynamic pruning result with its fixed core,
  plus core store/restore and invalid-selector checks.
- Two complete tree-likelihood tests:
  - all four selector values reproduce the corresponding fixed-model log
    likelihood to `1E-12`;
  - a proposed model change restores the original selector and likelihood
    exactly after rejection.
- Five prior/operator tests covering equal and unequal weights, zero weights,
  invalid weights, integer support, and the zero Hastings ratio.
- A 10,000-state smoke run with:

  ```text
  dr.app.beast.BeastMain -beagle_off -overwrite -seed 12345 \
      examples/release/flipflop/modelAveraging.xml
  ```

  The run completed in 0.296 seconds. All four selector values appeared in the
  log. The division-model operator made 9,900 proposals with acceptance
  probability 0.7921.

The existing `ant junit_flipflop` target retained its baseline result: eight
tests passed, while
`TestSubstitutionModelEmpiricalFrequencies` could not be discovered because
the source declares package `dr.evomodel.flipflop` but the Ant target requests
`test.dr.evomodel.flipflop`. This pre-existing unrelated test was not changed.

### Follow-up validation cleanup

After review, discarded calls to `getDivisionModelIndex` were replaced by the
explicit `validateDivisionModel` method. The discarded validation call was
removed from the division-model variable-change handler.

The selected pruning core is now cached. The handler refreshes that cache once
when `divisionModel` changes, validating the selector as it reads the new index,
then invalidates all nodes and fires the model-change event. Individual pruning
calls no longer read or validate the parameter. The current and stored selector
values follow the same MCMC store/restore behavior as the core's matrices and
partials. Because the parent tree likelihood skips core state storage when
`storePartials` is false, the model-averaging likelihood explicitly stores and
restores the selector cache in that mode.

The follow-up checks passed:

- `ant compile-all`.
- The three targeted model-averaging test classes: 10 tests passed, including
  selector rejection with both `storePartials` settings.
- `ant junit_flipflop`: all discoverable tests passed, with the same unrelated
  `TestSubstitutionModelEmpiricalFrequencies` class-discovery error described
  above.

## Deferred validation

Posterior model frequencies have not yet been compared with the existing
independent fixed-model simulation runs. That is the next validation phase.
