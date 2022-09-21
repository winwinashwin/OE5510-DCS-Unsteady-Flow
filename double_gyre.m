function [X,Y,phi,u_x,u_y] = double_gyre(x,y,t,A,e,w)
%DOUBLE_GYRE Compute stream function and velocity field for double gyre flow
[X,Y] = meshgrid(x,y);
if ~isa(t,'casadi.MX')
    t = reshape(t,1,1,[]);
end

a = e*sin(w*t);
b = 1 - 2*e*sin(w*t);
f = a.*X.^2 + b.*X;

phi = A*sin(pi*f).*sin(pi*Y);  % stream function

u_x = -pi*A*sin(pi*f).*cos(pi*Y);
u_y =  pi*A*cos(pi*f).*sin(pi*Y);
end
