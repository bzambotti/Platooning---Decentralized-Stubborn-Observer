# Distributed Stubborn Observer for Platooning

This project implements a Distributed Stubborn Observer (DSO) and compares its performance against a Centralized Stubborn Observer (CSO) within the context of connected vehicle platooning.

The code evaluates how a distributed system can reliably estimate the internal state of followers in a vehicle platoon using localized measurements, proving resilient against both spatial propagation of high-magnitude sporadic faults, and temporal accumulation of errors during periodic attacks.

## Scenarios
The simulation script `run_scenarios.m` executes the following test scenarios:
1. **Scenario 1 - Spatial Isolation:** Demonstrates that a localized measurement anomaly or high-magnitude fault (e.g., on a single vehicle sensor) does not propagate spatially to the rest of the interconnected decentralized observer nodes.
2. **Scenario 2 - Temporal Relaxation:** Demonstrates the ability of the decentralized model with finite memory buffers to quickly recover from periodic attack signals without accumulating significant state estimation error, contrasting with the Centralized observer.

## Prerequisites
To run this project, you need:
- **MATLAB** or **GNU Octave**.
- The **Control System Toolbox** installed (`pkg load control` for Octave).

## Project Structure
- `run_scenarios.m`: The primary simulation script.
- `src/`: Directory containing all supporting code and functions.
  - `model_platoon.m`: Derives the discrete-time state-space dynamics of the platoon.
  - `run_simulation.m`: Conducts the discrete-time simulation for a specified maneuver under disturbance.
  - `observer_stubborn_finite_memory.m`: Computes the centralized stubborn observer update.
  - `stubborn_finitememory.m`: Computes the distributed, decentralized finite memory stubborn observer update.
  - `generate_scenario_disturbance.m`: Utility functions for building attack profiles.
  - `generate_observer_gain.m`: Computes decoupled block LQR state observer gains.
  - `findAndTransformObserver.m`: Observability decomposition to isolate target node states.
  - `create_stubborn_matrices_robust.m`: Constructs mapping matrices for the stubborn memory buffer.

## Execution
Run the following script to instantiate the platoon dynamic model, build the distributed observer architectures, generate scenario specific disturbances, run the simulations, and render the evaluation plots.

### MATLAB
```matlab
run_scenarios
```

### GNU Octave
```bash
octave --eval "run_scenarios"
```
