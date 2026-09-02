# Hop-Challenge

Groundwork for a flight vehicle capable of launching and landing (a
hopper/VTVL vehicle), per the requirements in `Lander Challenge Rules.pdf`.

## Repo layout

- **`1D-throttle/`** — the only implemented project so far. A 1D powered-
  flight trajectory optimizer used to size the engine's thrust/throttle
  range and estimate propellant trade-offs for a hop. See below.
- `matlab-controls/` — not yet developed. Currently a copy of the
  Simulink dispersion model's scaffolding (`controls.slx`, `init_params.m`,
  reference tables), intended to become the vehicle's controls model.
- `NASA-CEA/` — not yet developed. A single script
  (`propellant_trade.py`) that uses `rocketcea` to tabulate Isp/C*/Tc/gamma
  across chamber pressure, mixture ratio, and area ratio for propellant
  selection.

## `1D-throttle/`

`main.py` sets up and solves a 1D vertical powered-flight trajectory
optimization for a hop — the ascent leg from the pad up to a target
altitude, arriving with zero velocity. The goal is to determine, for a
given engine's throttle range and Isp, the minimum-time thrust/throttle
profile and how much propellant it costs — inputs the docstring says will
feed into "necessary thrust and throttling capabilities" and "propellant
trade-offs" for the hopper, which then feed a full 6DOF controls program
(with TVC).

**Method:** direct (trapezoidal) collocation. State and control at `N + 1`
nodes, plus the free final time, are packed into one decision vector and
solved with `scipy.optimize.minimize` (SLSQP), warm-starting each stage of
a penalty continuation from the previous stage's solution.

- **State** `[h, v, m]` — altitude, velocity, vehicle mass (dry + remaining
  propellant).
- **Control** `[power, throttle]` — `power` gates the engine on/off,
  `throttle` sets commanded thrust between the engine's throttle floor
  (`T_MIN = THROTTLE * T_MAX`, currently 40% of max) and `T_MAX` when on.
- **Dynamics** — point-mass vertical equations of motion: thrust minus
  quadratic drag minus weight over mass for `v_dot`, and propellant
  depletion via `Isp` for `m_dot`.
- **Objective** — minimize final time, plus a penalty term
  (`power * (1 - power)`, ramped up over a continuation schedule of
  penalty weights `[0, 1, 10, 100, 1000]`) that pushes the `power` control
  toward bang-bang (fully on/off) rather than partial-throttle switching,
  since a real engine here is expected to ignite once rather than restart
  mid-flight.
- **Constraints** — trapezoidal dynamics defects between nodes, initial/
  final boundary conditions (start at the pad at rest with full propellant,
  arrive at the target altitude at rest), and a mass floor (never below dry
  mass).

**Run it:**

```bash
cd 1D-throttle
python main.py
```

Prints the SLSQP result at each penalty-continuation stage (final time,
number of still-fractional `power` nodes, objective, solver status), then
the converged final time and propellant used.

**Status:** `solve()` runs end-to-end and converges cleanly at every
continuation stage. `solve_reduced()` (meant to snap the optimized `power`
control to a hard single on/off ignition window and re-verify/re-solve
against it) and `dynamics_bounds()` are still unimplemented stubs — the
next step if this is picked back up.
