# Division-model averaging XML converter

## Summary

`simulate/convert_to_division_model_averaging.py` converts one fixed-model
PHYFUM XML file into an equivalent model-averaging configuration. It uses
targeted text transformations so comments and the surrounding hand-authored
layout remain intact, then parses the result to ensure it is well-formed XML.

The converter:

- replaces the fixed cenancestor likelihood and its references;
- adds the division-model parameter, operator, prior, and log entries;
- inserts `divisionModelAveraged` before the final extension of every
  `fileName` and `operatorAnalysis` value;
- initializes the identity model with a flat model prior and operator weight
  `0.5` by default; and
- removes power-posterior MLE estimation by default while preserving harmonic
  mean analysis.

Run `python3 simulate/convert_to_division_model_averaging.py --help` for all
options. A typical conversion is:

```text
python3 simulate/convert_to_division_model_averaging.py \
    -i fixed-model.xml \
    -o averaged-model.xml
```

Use bare `--MLE` or `--MLE 1` to retain an existing
`marginalLikelihoodEstimator`; `--MLE 0` and `--noMLE` remove it.

## Verification

The following checks passed on July 31, 2026:

- Nine focused Python unit and CLI tests covering defaults, short and long
  argument names, custom models and weights, MLE modes, output-name handling,
  invalid values, and already-converted input.
- Conversion and structural inspection of
  `sim_3_0.5_0.05_0.0625_100_5_1.3cells.xml`. Its original 750,000-state MCMC
  was not run.
- `ant compile-all`.
- A 20-state `-beagle_off` BEAST smoke run made from a temporary conversion of
  `examples/release/flipflop/sim1AbsoluteFission.xml`. The run parsed the new
  elements, selected the Java model-averaging likelihood core, logged the
  identity selector, and completed normally with seed `12345`.
