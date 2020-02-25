%Question 9
clc; clear;
%Givens
a = 20000;
e = 0.45;
Omega = deg2rad(59);
inc = deg2rad(27);
omega = deg2rad(94);
nu = deg2rad(58);

oe = [a; e; Omega; inc; omega; nu];
mu = 398600;

[rv,vv]  = oe2rv_Elosegui_Marcus(oe,mu);

fprintf('The position of the spacecraft: [%.4f\t%.4f\t%.4f] km\n',rv);
fprintf('The velocity of the spacecraft: [%.4f\t%.4f\t%.4f] km/s\n',vv);