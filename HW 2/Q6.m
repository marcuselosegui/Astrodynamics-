%Question 1-6
clc; clear;

mu = 1; 
rp = (1/20)* [-12; -20; 15];
a = -mu/(norm(rp))^3 * rp;
fprintf('The inertial acceleration of the spacecraft is [%.4f\t%.4f\t%.4f]\n',a(1),a(2),a(3))