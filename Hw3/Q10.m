%Question 10
clc; clear;
%Givens
a = 1.6;
e = 0.4;
Omega = deg2rad(287);
inc = deg2rad(46);
omega = deg2rad(28);
nu = deg2rad(139);

oe = [a; e; Omega; inc; omega; nu];
mu = 1;

[rv,vv]  = oe2rv_Elosegui_Marcus(oe,mu);

fprintf('The position of the spacecraft: [%.4f\t%.4f\t%.4f]\n',rv);
fprintf('The velocity of the spacecraft: [%.4f\t%.4f\t%.4f]\n',vv);