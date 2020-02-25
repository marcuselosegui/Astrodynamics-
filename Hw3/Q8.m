%Question 8
%Question 8
clc; clear;
%Givens
a = 19133.333;
e = 0.5;
Omega = deg2rad(30);
inc = deg2rad(45);
omega = deg2rad(45);
nu = deg2rad(0);

oe = [a; e; Omega; inc; omega; nu];
mu = 398600;

[rv,vv]  = oe2rv_Elosegui_Marcus(oe,mu);

fprintf('The position of the spacecraft: [%.4f\t%.4f\t%.4f] km\n',rv);
fprintf('The velocity of the spacecraft: [%.4f\t%.4f\t%.4f] km/s\n',vv);