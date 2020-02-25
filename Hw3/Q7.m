%Question 7
clc; clear;

%Givens
a = 15307.548;
e = 0.7;
Omega = deg2rad(194);
inc = deg2rad(39);
omega = deg2rad(85);
nu = deg2rad(48);

%Orbital Elements
oe = [a; e; Omega; inc; omega; nu];
mu = 398600;

[rv,vv]  = oe2rv_Elosegui_Marcus(oe,mu);

fprintf('The position of the spacecraft: [%.4f\t%.4f\t%.4f] km\n',rv);
fprintf('The velocity of the spacecraft: [%.4f\t%.4f\t%.4f] km/s\n',vv);