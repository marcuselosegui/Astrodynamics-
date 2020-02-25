%Question 6
clc; clear;

mu = 398600;
Re = 6371;

%orbit 1
h = 300;
r1 = Re+h;
inc1 = deg2rad(28.5);
e1 = 0;
omega1 = 0;
Omega1 = 0;
nu1 = 0;
oe1 = [r1;e1;Omega1;inc1;omega1;nu1];
[rv1,vv1] = oe2rv_Elosegui_Marcus(oe1,mu);
hv1 = cross(rv1,vv1);
h1 = norm(hv1);

%orbit 2
tauh = 24;
tau = 24*3600;
r2 = (mu*(tau/(2*pi))^2)^(1/3);
e2 = 0;
omega2 = 0;
Omega2 = Omega1;
inc2 = 0;
nu2 = 0;
oe2 = [r2;e2;Omega2;inc2;omega2;nu2];
[rv2,vv2] = oe2rv_Elosegui_Marcus(oe2,mu);
hv2 = cross(rv2,vv2);
h2 = norm(hv2);

a = (r1+r2)/2;
e = (r2-r1)/(r2+r1);

crank = dot(hv1,hv2)/(h1*h2);

v1m = sqrt(mu/r1);
v1p = sqrt(mu)*sqrt(2/r1 - 1/a);
v2m = sqrt(mu)*sqrt(2/r2 - 1/a);
v2p = sqrt(mu/r2);

%a
for f = 0:1
    dV1 = sqrt(v1m^2+v1p^2-2*v1m*v1p*cos(f*crank));
    dV2 = sqrt(v2m^2+v2p^2-2*v2m*v2p*cos((1-f)*crank));
    dV(f+1) = dV1+dV2;
end

f = [0:0.01:1];
vc1 = sqrt(mu/r1);

for i = 1:length(f)
    dVi1 = sqrt(v1m^2+v1p^2-2*v1m*v1p*cos(f(i)*crank));
    dVi2 = sqrt(v2m^2+v2p^2-2*v2m*v2p*cos((1-f(i))*crank));
    dVi(i) = dVi1+dVi2;
end
plot(f,dVi/vc1)

M = min(dVi);
n = find(dVi == M);
minf = f(n);

fprintf('a) The magnitude of the impulse with inclination change\n')
fprintf('  at apoapsis:    %.8f km/s\n',dV(1))
fprintf('b) The magnitude of the impulse with inclination change\n')
fprintf('  at initial LEO: %.8f km/s\n',dV(2))
fprintf('d) The value of f for smallest total impulse: %.8f\n',minf)