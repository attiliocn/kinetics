"""Gaussian perturbation mode (Run Mode 2) utilities."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from . import config


def perturb_energies(
    energies: pd.Series,
    sigma: float,
    rng: np.random.Generator,
    intermediates: Sequence[str] | None = None,
    transition_states: Sequence[str] | None = None,
) -> tuple[pd.Series, dict[str, float]]:
    """
    Apply one Gaussian perturbation sample to a mechanism-energy series.

    Sampling logic:
    1. A ~ N(0, sigma)
    2. For each transition-state label j, sample x_j ~ U(0, 1)
    3. For each transition-state label j, compute B_j = A * x_j
    4. intermediate labels receive +A
    5. each transition-state label j receives +B_j

    Args:
        energies: Baseline mechanism energies indexed by step label.
        sigma: Standard deviation for A.
        rng: NumPy random generator.
        intermediates: Labels treated as intermediates (defaults to config.INTERMEDIATES).
        transition_states: Labels treated as transition states
            (defaults to config.TRANSITION_STATES_LABELS).

    Returns:
        A tuple (perturbed_energies, sample_info), where sample_info contains
        scalar values plus per-transition-state x_j/B_j fields.
    """
    intermediates = list(config.INTERMEDIATES if intermediates is None else intermediates)
    transition_states = list(
        config.TRANSITION_STATES_LABELS if transition_states is None else transition_states
    )

    perturbed = energies.copy()

    a_sample = float(rng.normal(0.0, sigma))

    intermediate_idx = [label for label in intermediates if label in perturbed.index]
    transition_state_idx = [label for label in transition_states if label in perturbed.index]

    x_by_ts: dict[str, float] = {
        label: float(rng.uniform(0.0, 1.0)) for label in transition_state_idx
    }
    b_by_ts: dict[str, float] = {label: a_sample * x_val for label, x_val in x_by_ts.items()}

    if intermediate_idx:
        perturbed.loc[intermediate_idx] = perturbed.loc[intermediate_idx] + a_sample

    for label in transition_state_idx:
        perturbed.loc[label] = perturbed.loc[label] + b_by_ts[label]

    x_values = list(x_by_ts.values())
    b_values = list(b_by_ts.values())

    # Flatten per-label perturbation factors to CSV-friendly metadata columns.
    x_fields = {f"x_{label}": x_val for label, x_val in x_by_ts.items()}
    b_fields = {f"B_{label}": b_val for label, b_val in b_by_ts.items()}

    sample_info = {
        "A": a_sample,
        "x_mean": float(np.mean(x_values)) if x_values else float("nan"),
        "B_mean": float(np.mean(b_values)) if b_values else float("nan"),
        "n_intermediates_applied": float(len(intermediate_idx)),
        "n_transition_states_applied": float(len(transition_state_idx)),
        **x_fields,
        **b_fields,
    }
    return perturbed, sample_info


def aggregate_trajectories(
    trajectories: list[pd.DataFrame],
    time_column: str = "time",
    percentiles: Sequence[float] | None = None,
) -> dict[str, pd.DataFrame]:
    """
    Aggregate an ensemble of trajectory DataFrames.

    All trajectories must share the same time grid and species columns.

    Returns:
        {
          "mean": DataFrame,
          "std": DataFrame,
          "percentiles": DataFrame  # long format: time, species, percentile, value
        }
    """
    if not trajectories:
        raise ValueError("trajectories is empty")

    percentiles = (
        config.PERTURBATION_PERCENTILES_DEFAULT if percentiles is None else tuple(percentiles)
    )

    reference_cols = list(trajectories[0].columns)
    if time_column not in reference_cols:
        raise ValueError(f"Missing required time column '{time_column}' in trajectories")

    for idx, traj in enumerate(trajectories[1:], start=1):
        if list(traj.columns) != reference_cols:
            raise ValueError(
                f"Trajectory at index {idx} has different columns than index 0. "
                "Aggregation requires aligned trajectories."
            )
        if not np.array_equal(traj[time_column].to_numpy(), trajectories[0][time_column].to_numpy()):
            raise ValueError(
                f"Trajectory at index {idx} has a different time grid than index 0."
            )

    species_cols = [c for c in reference_cols if c != time_column]
    stacked = np.stack([traj[species_cols].to_numpy() for traj in trajectories], axis=0)

    mean_df = pd.DataFrame(stacked.mean(axis=0), columns=species_cols)
    std_df = pd.DataFrame(stacked.std(axis=0), columns=species_cols)
    mean_df.insert(0, time_column, trajectories[0][time_column].to_numpy())
    std_df.insert(0, time_column, trajectories[0][time_column].to_numpy())

    pct_blocks = []
    for p in percentiles:
        pct_values = np.percentile(stacked, p, axis=0)
        block = pd.DataFrame(pct_values, columns=species_cols)
        block.insert(0, time_column, trajectories[0][time_column].to_numpy())
        block = block.melt(id_vars=[time_column], var_name="species", value_name="value")
        block.insert(2, "percentile", p)
        pct_blocks.append(block)
    percentiles_df = pd.concat(pct_blocks, ignore_index=True)

    return {"mean": mean_df, "std": std_df, "percentiles": percentiles_df}


def run_reaction_pipeline_mode2(
    baseline_energies: pd.Series,
    run_single_pipeline: Callable[[pd.Series], pd.DataFrame],
    n_samples: int = config.PERTURBATION_N_SAMPLES_DEFAULT,
    sigma: float = config.PERTURBATION_SIGMA_DEFAULT,
    random_seed: int = config.PERTURBATION_SEED_DEFAULT,
    intermediates: Sequence[str] | None = None,
    transition_states: Sequence[str] | None = None,
    aggregate: bool = True,
    export_trajectories: bool = False,
    export_dir: str | Path | None = None,
    reaction_label: str | None = None,
) -> dict[str, Any]:
    """
    Execute Gaussian perturbation mode for one reaction.

    The provided `run_single_pipeline` callable must execute the full
    deterministic pipeline for one energy series and return a trajectory
    DataFrame with a `time` column and one column per species.

    This function always executes one unperturbed baseline run first, then
    executes `n_samples` perturbed runs. Therefore total runs are
    `n_samples + 1`.
    """
    if n_samples < 0:
        raise ValueError("n_samples must be >= 0")
    if sigma < 0:
        raise ValueError("sigma must be >= 0")
    if export_trajectories and (export_dir is None or reaction_label is None):
        raise ValueError(
            "When export_trajectories=True, both export_dir and reaction_label are required"
        )

    rng = np.random.default_rng(random_seed)

    trajectories: list[pd.DataFrame] = []
    sample_rows: list[dict[str, float]] = []
    perturbed_rows: list[dict[str, Any]] = []
    sample_ids: list[int] = []

    for run_idx in range(n_samples + 1):
        is_baseline = run_idx == 0
        if is_baseline:
            perturbed_energies = baseline_energies.copy()
            sample_info = {
                "A": 0.0,
                "x_mean": float("nan"),
                "B_mean": 0.0,
                "n_intermediates_applied": 0.0,
                "n_transition_states_applied": 0.0,
            }
            sample_id = -1
        else:
            perturbed_energies, sample_info = perturb_energies(
                baseline_energies,
                sigma=sigma,
                rng=rng,
                intermediates=intermediates,
                transition_states=transition_states,
            )
            sample_id = run_idx - 1

        trajectory = run_single_pipeline(perturbed_energies)
        if not isinstance(trajectory, pd.DataFrame):
            raise TypeError(
                "run_single_pipeline must return a pandas DataFrame trajectory"
            )
        if "time" not in trajectory.columns:
            raise ValueError("trajectory DataFrame must contain a 'time' column")

        trajectories.append(trajectory)
        sample_ids.append(sample_id)
        sample_rows.append(
            {
                "sample": float(sample_id),
                "is_baseline": float(is_baseline),
                **sample_info,
            }
        )
        perturbed_rows.append(
            {
                "sample": sample_id,
                "is_baseline": int(is_baseline),
                **{label: float(value) for label, value in perturbed_energies.items()},
            }
        )

    samples_df = pd.DataFrame(sample_rows)
    perturbed_energies_df = pd.DataFrame(perturbed_rows)
    result: dict[str, Any] = {
        "trajectories": trajectories,
        "samples": samples_df,
        "perturbed_energies": perturbed_energies_df,
    }

    if export_trajectories:
        out_dir = Path(export_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        trajectory_blocks: list[pd.DataFrame] = []
        for sample_id, traj in zip(sample_ids, trajectories, strict=True):
            block = traj.copy()
            block.insert(0, "sample", sample_id)
            block.insert(1, "is_baseline", int(sample_id == -1))
            block.insert(1, "reaction", reaction_label)
            trajectory_blocks.append(block)

        trajectories_long = pd.concat(trajectory_blocks, ignore_index=True)
        # CSV is very large for dense trajectory ensembles; use compressed binary storage.
        out_path = out_dir / f"{reaction_label}_mode2_trajectories.pkl.gz"
        trajectories_long.to_pickle(out_path, compression="gzip")
        result["trajectories_export_path"] = str(out_path)
        result["trajectories_export_format"] = "pickle_gzip"

        energies_out_path = out_dir / f"{reaction_label}_mode2_perturbed_energies.csv"
        perturbed_energies_df.to_csv(energies_out_path, index=False)
        result["perturbed_energies_export_path"] = str(energies_out_path)

    if aggregate and n_samples > 0:
        result["aggregate"] = aggregate_trajectories(trajectories[1:])

    return result
