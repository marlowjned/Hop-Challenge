'''
main.py
Throttle controls program for reaching target altitudes and landing.

Will be used for the hopper program to determine:
 - Necessary thrust and throttling capabilities
 - Propellant trade-offs

Which will then be used in a 6DOF controls progrram including TVC to ultimately 
determine the hopper system at a high level.

Author: Marlow Nedelchev
'''

import numpy as np
from scipy.optimize import minimize, NonlinearConstraint, Bounds

# === PARAMS === 
## VARIABLES
N = 20    # number of nodes

THROTTLE = 0.4 # Throttlability
T_MAX = 2000
T_MIN = T_MAX * THROTTLE

Isp = 200.0

# BCS         #  [h, v]
INITIAL_STATE =  [0, 0]
FINAL_STATE   = [50, 0]

DRY_MASS  = 25.0
PROP_MASS = 25.0
WET_MASS  = DRY_MASS + PROP_MASS

# GUESS
TF0 = 25.0
MASS_MARGIN = 5.0 # Prop mass left when reaching target

## TRUE CONSTANTS
NS = 3     # state dimension
NU = 2     # control dimension
NODE = NS + NU
NZ = (N + 1) * (NS + NU) + 1

g = 9.80665
CD = 0.5 # grab from ras
A_REF = np.pi * 0.1 ** 2
RHO_ATM = 1.225 # Assumed constant (dh < 100m)


# === DECISION VECTOR ===

def pack(X, U, Tf):
    z = np.zeros(NZ)
    z[:-1] = np.hstack([X, U]).reshape(-1)
    z[-1] = Tf
    return z

def unpack(z):
    blk = z[:-1].reshape(N + 1, NODE)
    X = blk[:, :NS]                  # (N+1, 3)
    U = blk[:, NS:NS + NU]           # (N+1, 2)
    tf = z[-1]
    return X, U, tf



# === DYNAMICS ===

def dynamics(state, controls):
    h, v, m = state[0], state[1], state[2]
    pow, thr = controls[0], controls[1]
    T = thrust(pow, thr)
    D = CD * A_REF * 0.5 * RHO_ATM * v * abs(v)
    hdot = v
    vdot = (T - D)/m - g
    mdot = -T/(Isp * g)
    return np.array([hdot, vdot, mdot])

# TODO: Run a version with single ignition
def thrust(power, throttle):
    return power * (T_MIN + throttle*(T_MAX - T_MIN))

# === OBJECTIVE FUNCTION ===

def objective_maker(penalty_prop):
    def objective(state):
        _, u, tf = unpack(state)
        power = u[:, 0]
        penalty = penalty_prop * np.sum(power * (1-power))
        return tf + penalty
    return objective


# === INITIAL GUESS ===

def init_guess():
    tf0 = TF0

    h_target = FINAL_STATE[0]
    v_avg = h_target / tf0
    h_lin = np.linspace(0.0, h_target, N + 1)
    v_lin = np.linspace(2*v_avg, 0.0, N + 1)

    m_lin = np.linspace(WET_MASS, DRY_MASS + MASS_MARGIN, N + 1)
    x0 = np.column_stack([h_lin, v_lin, m_lin])

    u0 = np.tile(np.array([1.0, 0.5]), (N + 1, 1))
    u0[int(0.5 * (N + 1)):, 0] = 0.0
    return pack(x0, u0, tf0)

def dynamics_bounds():
    pass


# === BOUNDS ===
# variable by variable constraints

def bounds():
    lb = np.full(NZ, -np.inf)
    ub = np.full(NZ,  np.inf)
    for k in range(N + 1):
        base = k * NODE
        # h >= 0
        lb[base + 0] = 0.0
        # v unbounded
        # m: dry <= m <= m0
        lb[base + 2] = DRY_MASS
        ub[base + 2] = WET_MASS
        # s in [0,1], u in [0,1]
        lb[base + NS + 0] = 0.0
        ub[base + NS + 0] = 1.0
        lb[base + NS + 1] = 0.0
        ub[base + NS + 1] = 1.0
    # t_f bounds
    lb[-1] = 1.0
    ub[-1] = 200.0
    return Bounds(lb, ub)

# === CONSTRAINTS ===
# constraints which are functions of potentially many variables

def dynamics_constraint(z):
    ''' xk+1 - xk - h/2 * (fk+1 - fk) = 0 '''
    X, U, tf = unpack(z)
    h_step = tf / N
    f = np.array([dynamics(X[k], U[k]) for k in range(N + 1)])
    d = X[1:] - X[:-1] - 0.5 * h_step * (f[:-1] + f[1:])
    return d.reshape(-1)

def boundary_constraint(z):
    ''' initial and final states '''
    X, U, tf = unpack(z)
    bc0 = np.array([X[0, 0] - INITIAL_STATE[0],
                    X[0, 1] - INITIAL_STATE[1],
                    X[0, 2] - WET_MASS])       
    bcf = np.array([X[-1, 0] - FINAL_STATE[0], 
                    X[-1, 1] - FINAL_STATE[1]])
    return np.concatenate([bc0, bcf])

def mass_constraint(z):
    ''' m >= m_dry '''
    X, _, _ = unpack(z)
    return X[:, 2] - DRY_MASS


# === SOLVER === 

def solve():
    guess = init_guess()
    bnds = bounds()
    cons = [
        {"type": "eq",   "fun": dynamics_constraint},
        {"type": "eq",   "fun": boundary_constraint},
        {"type": "ineq", "fun": mass_constraint},
    ]
 
    pow_penalty_vals = [0.0, 1.0, 10.0, 100.0, 1000.0]
    z = guess
    for penalty_val in pow_penalty_vals:
        res = minimize(
            objective_maker(penalty_val),
            z,
            method="SLSQP",
            bounds=bnds,
            constraints=cons,
            options={"maxiter": 300, "ftol": 1e-7, "disp": False},
        )
        z = res.x
        X, U, tf = unpack(z)
        s = U[:, 0]
        frac = np.sum((s > 1e-3) & (s < 1 - 1e-3))
        print(f"  rho={penalty_val:>7.1f}  tf={tf:7.3f}s  "
              f"fractional s nodes={frac:2d}  obj={res.fun:.4f}  "
              f"status={res.message}")
    return z

def solve_reduced(z):
    X, U, tf = unpack(z)
    pow_int = (U[:, 0] > 0.5).astype(float)


    pass


if __name__ == "__main__":
    z_final = solve()
    X, U, tf = unpack(z_final)
    prop_used = WET_MASS - X[-1, 2]
    print(f"\nFinal: tf={tf:.3f}s  propellant used={prop_used:.3f} kg")

