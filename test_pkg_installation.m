clc; close all; clear all;

%% Test CVX Installation
%      Simple convex optimization problem from compressed sensing context

n = 100;
p = 20;
Theta = randn(p,n);  
y = randn(p,1);

cvx_begin quiet;
    variable s_L1(n); 
    minimize( norm(s_L1,1) ); 
    subject to 
        Theta*s_L1 == y;
cvx_end;

fprintf("[+] CVX installation check complete\n");

%% Test CasADi Installation
%      Optimization problem HS071 (Hock & Schittkowski Problem #71)

import casadi.*

opti = Opti();

x1 = opti.variable();
x2 = opti.variable();
x3 = opti.variable();
x4 = opti.variable();

opti.minimize( x1*x4*(x1 + x2 + x3) + x3 );
opti.subject_to( x1*x2*x3*x4 >= 25 );
opti.subject_to( x1^2 + x2^2 + x3^2 + x4^2 == 40);
opti.subject_to( {
    1 <= x1 <= 5,
    1 <= x2 <= 5,
    1 <= x3 <= 5,
    1 <= x4 <= 5,
    } );

opti.set_initial(x1,1);
opti.set_initial(x2,5);
opti.set_initial(x3,5);
opti.set_initial(x4,1);

opts = struct;
opts.ipopt.print_level = 0;
opti.solver('ipopt',opts);

sol = opti.solve();

assert(abs(sol.value(x1)-1.00000000) < 1e-6);
assert(abs(sol.value(x2)-4.74299963) < 1e-6);
assert(abs(sol.value(x3)-3.82114998) < 1e-6);
assert(abs(sol.value(x4)-1.37940829) < 1e-6);

fprintf("[+] CasADi installation check complete\n");
