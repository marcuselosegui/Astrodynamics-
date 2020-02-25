%Question 1-10
clc; clear;

%givens
mu = 3.986004418e14;
me = -2e8 * 0.092903; %mechanical energy converted to m^2/s^2
e = 0.2;

a = -mu/(2*me);
p = a*(1-e^2);
h = sqrt(mu*p);
rp = a*(1-e);
ra = a*(1+e);

fprintf('(a)The magnitude of the specific angular momentum: %.4f m/s\n',h)
fprintf('(b)The semi-latus rectum: %.4f m\n',p)
fprintf('(c)The semi-major axis: %.4f m\n',a)
fprintf('(d)The periapsis radius: %.4f m\n',rp)
fprintf('(e)The apoapsis radius: %.4f m\n',ra)
