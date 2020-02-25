%Question 1-12
clc; clear;

%givens
mu = 398600;
re = 6378.145;
h = 4000;
v = 0.8; %changed to km/s
r = h+re;
gamma = 0;

me = (v/2) - mu/r; %mechanical energy
h = r*v*cos(gamma);
p = h^2/mu;

a = 1/((2/r)-(v^2/mu));
e = sqrt(1-(p/a));
rp = a*(1-e);
ra = a*(1+e);

fprintf('(a)The specific angular energy: %.4f km^2/s^2\n',me)
fprintf('(b)The magnitude of the specific angular momentum: %.4f km/s\n',h)
fprintf('(c)The semi-latus rectum: %.4f km\n',p)
fprintf('(d)The periapsis radius: %.4f km\n',rp)
fprintf('(e)The apoapsis radius: %.4f km\n',ra)