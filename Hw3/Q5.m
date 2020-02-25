%Question 5
clc; clear;

%givens
mu = 1;
rv = [0.7; 0.6; 0.3];
vv = [-0.8; 0.8; 0];

oe = rv2oe_Elosegui_Marcus(rv,vv,mu);

tau = 2*pi*sqrt(oe(1)^3/mu); %time period
p = oe(1)*(1-oe(2)^2); %semi-latus rectum
hv = cross(rv,vv); %The specific angular momentum
h = sqrt(p*mu); %Magnitude of the specific angular momentum
E = -mu/(2*oe(1)); %Specific mechanical energy

fprintf('(a) The six orbital elements for the orbit of the asteroid are:\n')
fprintf('-semi-major axis: %.4f AU\n',oe(1))
fprintf('-eccentricity: %.4f \n',oe(2))
fprintf('-longitude of the ascending node: %.4f rad\n',oe(3))
fprintf('-inclination: %.4f rad\n',oe(4))
fprintf('-argument of the periapsis: %.4f rad\n',oe(5))
fprintf('-true anomaly: %.4f rad\n',oe(6))
fprintf('(b) The orbital period: %.4f TU\n',tau)
fprintf('(c) The semi-latus rectum: %.4f AU\n', p)
fprintf('(d) The specific angular momentum: [%.4f\t%.4f\t%.4f] AU^2/TU\n',hv)
fprintf('-magnitude of the specific angular momentum: %.4f AU^2/TU\n',h)
fprintf('(e) The specific mechanical energy: %.4f AU^2/TU^2\n',E)