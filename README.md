# Dynamic Compressed Sensing of Periodic Unsteady Flow Fields

Given a periodic unsteady flow field and an environment parameter to observe (velocity, temperature, salinity etc), this work investigates a methodology that involves the following,
- Identify an optimal set of sensing waypoints via compressed sensing that minimizes the reconstruction error on proper orthogonal decomposition modes
- Optimize for an efficient sampling trajectory traversing these waypoints while minimizing actuation energy of the sensing robot and overall sensing duration via nonlinear optimization.

The implementation is carried out for an unsteady double-gyre flow field.

Heavily inspired by the paper [*"Dynamic Compressed Sensing of Unsteady Flows with a Mobile Robot"*](https://arxiv.org/abs/2110.08658)

## External Dependencies
- [CVX](http://cvxr.com/cvx/download/)
- [CasADi v3.5.5](https://web.casadi.org/get/)
