%Question 1-11
clc; clear;

%givens
mu = 398600; %mu in km taken from Q14
hp = 370;
e = 0.1;
re = 6378.145;

rp = hp+re;
a = rp/(1-e);
ra = a*(1+e);
ha = ra-re;

me = -mu/(2*a); %Mechanical Energy
p = a*(1-e^2);
h = sqrt(mu*p);

fprintf('(a)The apoapsis altitude: %.4f km\n',ha)
fprintf('(b)The specific mechanical energy: %.4f km^2/s^2\n',me)
fprintf('(c)The magnitude of the specific angular momentum: %.4f km/s\n',h)
fprintf('(d)The semi-latus rectum: %.4f km\n',p)