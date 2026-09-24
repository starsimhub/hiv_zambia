# Calibration experiments

Two-level structure: **epochs** and **experiments**.

## Epoch

A directory `experiments/NN_slug/` groups experiments that share the same
structural model setup. A structural change — a new data preprocessing step,
a new module, a change to the observation model, a change to the sim
duration/timestep — starts a new epoch. Experiments within an epoch are
directly comparable to each other; experiments across epochs are not.

Each epoch folder has a `README.md` describing what makes it distinct from
the previous epoch, so a future reader knows why the numbers can't be
compared across the boundary.

## Experiment

A directory `experiments/NN_epoch/NN_slug/` records a single calibration run
inside an epoch. Between experiments in the same epoch, only the
`calib_pars` dict changes (or the number of trials / workers). The model,
data files, and shared driver code are the same.

Each experiment folder contains:

- `run.py` — small driver: docstring describing this experiment, the
  `calib_pars` dict, and a call to `run_and_save` from the shared
  `run_hiv_calibration.py` at repo root. Runnable from the repo root as
  `python experiments/NN_epoch/NN_slug/run.py`.
- `SUMMARY.md` — date, commit hash, `calib_pars` verbatim, fit table,
  observations, next-step candidate. Embeds figures inline (`![](figures/…)`).
- `figures/` — snapshot of the plots produced from this experiment's run.

Raw calibration outputs (`raw_results/*.obj`) are gitignored bulk; the
reproducibility contract is `commit hash + run.py`. Between commits, the
shared driver (`run_hiv_calibration.py`) may evolve, and per-experiment
`run.py` files pin their own `calib_pars` — checking out an old commit and
running the folder's `run.py` reproduces that experiment.

## Deviation from calib-plugin default validator

The `calib:experiment-close` validator expects a `config.yaml` alongside
`run.py`. We do not use `config.yaml`: `calib_pars` lives verbatim in
`run.py` as a Python dict, which is the actual input to the calibration.
A yaml would be a redundant restatement. The other validator checks
(SUMMARY exists and is non-trivial, referenced figures resolve, run script
exists) are honored.
