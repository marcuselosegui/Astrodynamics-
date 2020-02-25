%Code for Question 1-3
clc;clear;

Re = 6378.145; %Earths Radius
rp = Re + 500; %Periapsis altitude
ra = Re + 800; %Apoapsis altitude
mu = 398600;

%Semi-Major Axis
a = (rp+ra)/2;

%Eccentrcity 
e = (ra - rp)/(ra + rp);

%semi-latus rectum
p = a*(1-e^2);

%magnitude of the specific angular momentum
h = sqrt(p*mu);

%speed of the spacecraft at periapsis and apoapsis
vp = sqrt(mu)*sqrt((2/rp) - (1/a));
va = sqrt(mu)*sqrt((2/ra) - (1/a));

fprintf('(a)The semi-major axis: %f km\n', a)
fprintf('(b)The eccentricity:  %f\n', e)
fprintf('(c)The semi-latus rectum %f km\n',p)
fprintf('(d)Magnitude of the specific angular momentum: %f km/s\n',h)
fprintf('(e)Speed of the spacecraft at periapsis: %f km/s\n',vp)
fprintf('(e)Speed of the spacecraft at apoapsis: %f km/s\n',va)