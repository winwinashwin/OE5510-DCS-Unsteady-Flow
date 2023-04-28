# Dynamic Compressed Sensing of Periodic Unsteady Flow Fields

Given a periodic unsteady flow field and an environment parameter to observe (velocity, temperature, salinity etc), this work investigates a methodology that involves the following,
- Identify an optimal set of sensing waypoints via compressed sensing that minimizes the reconstruction error on proper orthogonal decomposition modes
- Optimize for an efficient sampling trajectory traversing these waypoints while minimizing actuation energy of the sensing robot and overall sensing duration via nonlinear optimization.

The implementation is carried out for an unsteady double-gyre flow field.

Heavily inspired by the paper [*"Dynamic Compressed Sensing of Unsteady Flows with a Mobile Robot"*](https://arxiv.org/abs/2110.08658)

## External Dependencies
- [CVX](http://cvxr.com/cvx/download/)
- [CasADi v3.5.5](https://web.casadi.org/get/)

### Installation

#### CVX
- Download and extract the _Redistributable, all platform version_ of CVX v2.2 to a folder from the [official download page](http://cvxr.com/cvx/download/)
- Navigate to the CVX folder from MATLAB command line and run `cvx_setup` command

#### CasADi
- Download and extract CasADi v3.5.5 corresponding to your OS and MATLAB version from the [official download page](https://web.casadi.org/get/)
- From MATLAB command line, add CasADi to MATLAB path by running the command
```MATLAB
addpath('<your-casadi-installation-path>')
```

Please run the test script `src/test_pkg_installation.m` to make sure the packages have been installed and setup correctly.
The script should run and exit without any error upon successful setup.
