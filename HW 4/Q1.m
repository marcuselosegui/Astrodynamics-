%Question 1
clc;clear;
mu = 398600;
tspan = [21.02*60,1913.38*60];
z = [-664.699;8112.75;4479.81;-0.87036;-0.068046;-8.290459];
options = odeset('Reltol',1e-6);  
[t,z] = ode113(@TwoBodyODE, tspan, z, options);
plot3(z(:,1),z(:,2),z(:,3))
fprintf('final PCI positon: [%.8f\t%.8f\t%.8f] km\n',z(end,1:3))
fprintf('final PCI positon: [%.8f\t%.8f\t%.8f] km/s\n',z(end,4:6))