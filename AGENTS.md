# Microkinetics Pipeline Implementation Plan

## Project Overview
Refactor the microkinetics notebook into modular Python packages with a clean notebook orchestration layer. The system supports two run modes:
- **Run Mode 1 (Baseline)**: Standard simulation with raw barrier energies
- **Run Mode 2 (Gaussian Perturbation)**: Uncertainty quantification via repeated sampling with randomly perturbed energies

---

## Phase 0: Project Setup & Design Decisions
**Status**: ✅ Completed

### 0.1 - Classify Species Types
**Description**: Determine which species in the reaction network are intermediates vs. transition states.

**Current assumption from original outline**:
- Intermediates: species where energy is directly perturbed by Gaussian sample A
- Transition states: species where energy is perturbed by B = A × x (where x ~ N(0, 1))

**Questions to resolve**:
- How is this classification encoded in the dataset (data structure, naming convention, separate metadata)?
- Are all 5 entities ['bispyr', 'roh', 'me2pyr', 'map', 'bis'] classified, or only a subset?
- Should classification be configurable (e.g., in `config.py`)?

**Tasks**:
- [x] Document how to distinguish intermediates from transition states in the current dataset
- [x] Create a mapping/classifier in `src/config.py`
- [x] Determine if classification is per-reaction or global

**Notes**:
- intermediates and Transition States will be identified by their label in the mechanism pd.Series, using a variable TRANSITION_STATES_LABELS:list and INTERMEDIATES:list
- We should have a cell at the top to set all relevant parameters like filepath and labels.
- All entries should be classified either as TS or intermediate

---

### 0.2 - Define "Mechanism" Representation
**Description**: Clarify what the "mechanism" output should be (energy profile, reaction pathway, or other).

**Current assumption**: 
- A sequence or structure representing the energy landscape for the reaction

**Questions to resolve**:
- Is "mechanism" a pd.Series indexed by reaction step?
- Does it include information about barriers, intermediates, or both?
- Should it be derived from the energies or is it a fixed property of each reaction?

**Tasks**:
- [x] Confirm mechanism representation (pd.Series, dict, custom object)
- [x] Determine if mechanism is computed from energies or looked up
- [x] Document the exact data flow: energies → mechanism → activation energies

**Notes**:
- A mechanism is a pd.Series with labels=reaction steps (reactants, int1 etc) and values are relative energies (reactants = 0)
- Only include the relative energies, derived from the energies. A function will map to a particular subset=series from the energies dataframe

---

### 0.3 - Finalize Perturbation Logic
**Description**: Lock down the exact mathematical procedure for Run Mode 2.

**Current specification from outline**:
```
For each iteration n in [1, n_samples]:
  1. Sample A ~ N(0, sigma)
  2. Sample x ~ U(0, 1), compute B = A * x
  3. For each intermediate: perturbed_energy = energy + A
  4. For each transition state: perturbed_energy = energy + B
  5. Run full pipeline (energies → mechanism → ... → simulation)
  6. Store trajectory and any metrics
```

**Questions to resolve**:
- Should sigma be a single global scalar or per-species?
- Should the random seed be fixed for reproducibility, or configurable?
- What output should be collected from each iteration (full trajectory, summary stats, both)?
- Should we save all n trajectories or just aggregated results (mean, std, percentiles)?

**Tasks**:
- [x] Confirm perturbation sampling logic (especially how A and B are used)
- [x] Define what "each energy of an intermediate/transition state" means in this reaction network
- [x] Decide on output format (trajectories, aggregated stats, or both)
- [x] Specify random seed handling

**Notes**:
- Sigma is global, defaults to 2
- Random seed fixed
- For each iteration we should have a dataframe of time and several columns for each label containing its concentration at a particular time. This is going to be massive so some kind of optimization might be needed instead of a list of dataframes
- We save all trajectories and analyze the data using other cells specific for that

---

## Phase 1: Core Modules (No Perturbation Yet)
**Status**: ✅ Completed

These modules implement Run Mode 1 (baseline). Perturbation logic is added in Phase 2.

### 1.1 - `src/config.py`
**Description**: Central configuration and constants for the reaction system.

**Contents**:
- `ENTITIES`: List of all species in the network (e.g., `['bispyr', 'roh', 'me2pyr', 'map', 'bis']`)
- `INTERMEDIATE_SPECIES`: Subset of ENTITIES that are intermediates (to be confirmed in Phase 0.1)
- `TRANSITION_STATE_SPECIES`: Subset of ENTITIES that are transition states (to be confirmed in Phase 0.1)
- `INITIAL_CONCENTRATIONS`: Default starting concentrations (dict by species name)
- `ODE_SOLVER_OPTIONS`: Dict with method, rtol, atol, etc.

**Tasks**:
- [x] Create file and import it successfully in main.ipynb
- [x] Define ENTITIES and verify against actual data
- [x] Document species classification (intermediates vs. transition states)
- [x] Add version/notes explaining any hardcoded assumptions

**Status**: ✅ Completed

---

### 1.2 - `src/kinetics.py`
**Description**: Kinetic theory calculations (barrier → rate constant).

**Functions** (move from current notebook):
- `calculate_rate(reaction_barrier: float, temperature: float) -> float`
  - Input: barrier in kcal/mol, temperature in Celsius
  - Output: rate constant (s⁻¹)
  - Physics: Arrhenius-like calculation using fundamental constants

**Tests to add**:
- [ ] Verify units (kcal/mol input → s⁻¹ output at reference T)
- [ ] Check edge cases (T=0, very high barriers, etc.)

**Tasks**:
- [x] Create file with `calculate_rate()` function
- [x] Add docstring with unit details
- [x] Test against known barrier/rate pairs
- [x] Import and verify in notebook

**Status**: ✅ Completed

---

### 1.3 - `src/ode_system.py`
**Description**: ODE system definition for the reaction mechanism.

**Functions**:
- `kinetic_equations(t: float, c: list, *rates) -> list`
  - Input: time, concentration vector [bispyr, roh, me2pyr, map, bis], rate constants
  - Output: time derivatives dc/dt for each species
  - Math: Hard-coded mass-action kinetics for the specific reaction network

**Tasks**:
- [x] Move `kinetic_equations()` from notebook
- [x] Verify entity order matches `config.ENTITIES`
- [x] Add helper to convert between list ↔ dict representations
- [x] Document the reaction scheme (which steps, which rates)
- [x] Test with known rate constants + initial conditions

**Status**: ✅ Completed

---

### 1.4 - `src/simulation.py`
**Description**: Numerical integration and trajectory utilities.

**Functions**:
- `run_simulation(rates: dict, concentrations: dict, times: list) -> OdeResult`
  - Wraps `scipy.integrate.solve_ivp()` with configured ODE system
  - Uses options from `config.ODE_SOLVER_OPTIONS`
- `exponential_trajectory_points(x_end: float, num_points: int = 100, exponent: float = 3) -> np.ndarray`
  - Generate time points clustered at the start

**Tasks**:
- [x] Move `run_simulation()` and `exponential_trajectory_points()` from notebook
- [x] Import `kinetic_equations` from `src/ode_system.py`
- [x] Test simulations produce expected concentration profiles
- [x] Verify time-point distribution for various exponents

**Status**: ✅ Completed

---

### 1.5 - `src/data.py`
**Description**: Data loading, parsing, filtering.

**Functions**:
- `load_dataset(filepath: str, sheet_name: str = 'barriers') -> pd.DataFrame`
  - Load Excel file with barrier energies
- `parse_reaction_id(reaction_label: str) -> dict`
  - Extract imido and roh type from label like "im1-w1m"
- `select_reactions(df: pd.DataFrame, filters: dict[str, Any]) -> pd.DataFrame`
  - Apply boolean filters (e.g., `{'imido': 'im1', 'roh': ['w1m', 'd1m', ...]}`)
- `compute_rate_constants(barriers: pd.DataFrame, temperature: float = 25) -> pd.DataFrame`
  - Apply `calculate_rate()` to all barriers, rename columns ('b1' → 'k1d', etc.)

**Tasks**:
- [x] Create file and move dataset loading logic from notebook
- [x] Test with actual `calculated-exchange.xlsx`
- [x] Verify reaction ID parsing handles all expected formats
- [x] Verify filtering returns correct subset

**Status**: ✅ Completed

---

### 1.6 - `src/plotting.py`
**Description**: Visualization utilities.

**Functions**:
- `plot_kinetics_roh(times: np.ndarray, concentrations: dict, reaction_label: str, colors=None, **kwargs) -> None`
  - Create and save concentration vs. time plot
  - **Key change**: colors should be passed in, not global

**Tasks**:
- [x] Move `plot_kinetics_roh()` from notebook
- [x] Make colors parameterized (not global)
- [x] Add optional toggle for ylim/xlim settings
- [x] Test plot generation and file saving

**Status**: ✅ Completed

---

### 1.7 - Main Notebook Orchestration (Run Mode 1)
**Description**: Rewrite the main notebook to be a clean 4-cell orchestration script.

**Cell structure**:

**Cell 1: Setup & Imports**
- Import all modules from `src/`
- Setup matplotlib/seaborn
- Define color palette and display config
- Status: ✅ Ready (modules exist)

**Cell 2: Load & Filter Data**
- Load dataset using `src.data.load_dataset()`
- Parse reaction IDs and apply filters
- Compute rate constants
- Status: ✅ Ready (functions available)

**Cell 3: Run Mode 1 (Baseline)**
- Loop through selected reactions
- For each reaction, call high-level pipeline function
- Collect results (trajectories, conversions, yields)
- Status: ✅ Ready (core functions available)

**Cell 4: Export & Report**
- Save results to CSV
- Generate summary statistics
- Optional: create plots
- Status: ✅ Ready (plotting available)

---

## Phase 2: Gaussian Perturbation Mode
**Status**: ✅ Completed

Extends Phase 1 modules to support Run Mode 2 (randomized energy perturbations).

### 2.1 - Design: Perturbation Logic
**Description**: Plan how perturbations fit into the pipeline.

**Key decisions needed** (from Phase 0):
- How to identify intermediates vs. transition states
- How to apply A and B to energies before recomputing rate constants

**Proposed structure**:
- New function `src/perturbation.py` with:
  - `perturb_energies(energies: dict, sigma: float, rng: np.random.Generator) -> dict`
    - Sample A, x, compute B
    - Apply to appropriate species
    - Return perturbed energies

**Tasks**:
- [x] Finalize species classification (Phase 0.1)
- [x] Define perturbation algorithm precisely
- [x] Create `src/perturbation.py` with perturbation logic
- [x] Unit test perturbation on simple cases

**Status**: ✅ Completed

---

### 2.2 - Perturbation Pipeline Integration
**Description**: Wire perturbation into the simulation loop.

**Changes**:
- Wrapper function: `run_perturbation_mode(baseline_energies, n_samples, sigma, rng_seed)`
  - Calls `perturb_energies()` in a loop
  - For each sample, runs full baseline pipeline on perturbed energies
  - Collects all trajectories

**Tasks**:
- [x] Create `run_reaction_pipeline_mode2()` function that wraps Mode 1 with perturbation
- [x] Test on single reaction, verify n_samples iterations run
- [x] Handle random seed for reproducibility

**Status**: ✅ Completed

---

### 2.3 - Result Aggregation for Mode 2
**Description**: Summarize perturbed trajectory ensemble.

**New functions** in `src/` (TBD which module):
- `aggregate_trajectories(trajectories: list[pd.DataFrame]) -> dict`
  - Compute mean, std, percentiles (5, 25, 50, 75, 95) for each species
  - Return as DataFrame or nested dict

**Tasks**:
- [x] Define aggregation metrics (mean, std, percentiles, etc.)
- [x] Implement aggregation function
- [x] Test on synthetic trajectory ensemble

**Status**: ✅ Completed

---

### 2.4 - Notebook Cell: Run Mode 2
**Description**: Add notebook cell for Gaussian perturbation execution.

**Cell structure**:
- Input: n_samples (default 1000), sigma (default TBD), random seed
- Loop through selected reactions
- For each: call Mode 2 pipeline, collect ensemble
- Aggregate results, display alongside Mode 1
- Status: ✅ Ready (utilities implemented)

---

## Phase 3: Comparison & Analysis
**Status**: ⏸ Not Planned (Out of Scope)

### 3.1 - Side-by-side Mode 1 vs Mode 2
**Description**: Compare baseline and perturbed results.

**Tasks**:
- [ ] Create comparison DataFrame (baseline conversion vs. Mode 2 mean ± std)
- [ ] Visualize baseline vs. ensemble trajectories on same plot
- [ ] Export comparison table

**Status**: ⏸ Not Planned

---

### 3.2 - Sensitivity Analysis (Optional)
**Description**: Vary sigma and compare outputs.

**Tasks**:
- [ ] Run Mode 2 for multiple sigma values
- [ ] Plot output (e.g., conversion) vs. sigma
- [ ] Identify which species are most sensitive

**Status**: ⏸ Not Planned

---

## Implementation Checklist

**Quick Reference** (copy-paste for progress tracking):
```
Phase 0 (Planning):
  [x] 0.1: Species classification
  [x] 0.2: Mechanism representation
  [x] 0.3: Perturbation logic

Phase 1 (Core):
  [x] 1.1: config.py
  [x] 1.2: kinetics.py
  [x] 1.3: ode_system.py
  [x] 1.4: simulation.py
  [x] 1.5: data.py
  [x] 1.6: plotting.py
  [x] 1.7: Notebook (Run Mode 1)

Phase 2 (Perturbation):
  [x] 2.1: Perturbation logic design
  [x] 2.2: Pipeline integration
  [x] 2.3: Result aggregation
  [x] 2.4: Notebook (Run Mode 2)

Phase 3 (Analysis):
  [-] 3.1: Comparison tools (out of scope)
  [-] 3.2: Sensitivity (optional, out of scope)
```

---

## Notes & Decisions Log

**Date: Apr 7, 2026**
- Initial refactoring plan created
- Identified Phase 0 open questions (species classification, mechanism representation, perturbation details)
- Proposed modular structure with 6 core modules + notebook orchestration

**Date: Apr 7, 2026 (Phase 0 Complete)**
- Species classification: Global, binary classification (intermediate vs. transition state) via INTERMEDIATES and TRANSITION_STATES_LABELS lists in config
- Mechanism representation: pd.Series with reaction steps as index, relative energies as values (reactants = 0)
- Perturbation logic: Global sigma (default=2), fixed random seed, save all trajectories (each iteration returns full trajectory DataFrame)
- Phase 0 design decisions locked in; Phase 1 can now proceed with module implementation

**Date: Apr 7, 2026 (Phase 1 Complete)**
- Created src/ directory with 6 core modules: config.py, kinetics.py, ode_system.py, simulation.py, data.py, plotting.py
- Moved all chemistry functions from notebook to modules with full docstrings
- Added helper functions for concentration list/dict conversion
- config.py now defines all constants and species classification in one place
- All Phase 1 modules are tested and ready for notebook integration
- Phase 2 (Gaussian perturbation) can now be designed and implemented

**Date: Apr 7, 2026 (Phase 2 Complete)**
- Added `src/perturbation.py` with perturbation sampling logic `A ~ N(0, sigma)`, `x ~ U(0,1)`, `B = A*x`
- Implemented `perturb_energies()` using global label classification from config
- Implemented `run_reaction_pipeline_mode2()` wrapper to run full baseline pipeline over `n_samples`
- Implemented `aggregate_trajectories()` with mean/std/percentiles (5,25,50,75,95)
- Added a high-level Mode 2 example cell in the notebook for orchestration

**Date: Apr 7, 2026 (Phase 2 Update)**
- Added export of all sampled trajectories to CSV for each reaction
- Added tracking/export of per-sample perturbed mechanism energies for cross-checking (`*_mode2_perturbed_energies.csv`)
- Kept sample metadata export with A, x, and B for auditability

**Date: Apr 7, 2026 (Scope Decision)**
- Phase 3 (comparison/sensitivity) will not be implemented in this workstream
