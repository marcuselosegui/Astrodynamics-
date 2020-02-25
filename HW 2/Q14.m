%Question 1-14
clc;clear;

%givens
r = 403000;
mu = 398600;
v = 2.25;
nu = deg2rad(151);
re = 6378.145;

a = 1/(2/r - v^2/mu); %negative because hyperbolic orbit
x = [a r*cos(nu) (r-a)];
ecc = roots(x); %e must be >1

e = 0;

%used to sort through the roots to find the e>1
for i = 1:2
if ecc(i)>1
    e=ecc(i);
    i=i+1;
end
end

rp = a*(1-e);
hp = rp-re; %periapsis altitude

vp = sqrt(mu)*sqrt((2/rp) - (1/a));

fprintf('(a)The eccentricity: %.4f\n',e)
fprintf('(b)The periapsis altitdue: %.4f km\n',hp)
fprintf('(c)The periapsis speed: %.4f km/s\n',vp)

