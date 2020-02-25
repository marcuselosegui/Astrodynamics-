%Question 4
clc; clear;

mu = 398600;
rv = [2721.965;3522.863;5267.244];
vv = [9.572396;-0.474701;-2.725664];
t0 = 3.93*60;
tf  = 1771.58*60;

[r,v,E,nu] = propagateKepler_Elosegui_Marcus(rv,vv,t0,tf,mu);

%part e
n = 300; %increased from 100 to get a smoother curve
tinterval = transpose(linspace(t0,tf,n)); %creates 100 subintervals

RR = zeros(n,3);
VV = zeros(n,3);
for i = 1:n
    [RR(i,:),VV(i,:)] = propagateKepler_Elosegui_Marcus(rv,vv,t0,tinterval(i),mu); 
    %position and velocity at every time in tinterval
end
z = [transpose(rv) transpose(vv)]; %stores original position and 
%velocity in one row matrix

options = odeset('Reltol',1e-6);
[t,zv] = ode113(@TwoBodyODE, tinterval, z, options);


plot3(RR(:,1),RR(:,2),RR(:,3)) %plots from Kepler
hold on
plot3(zv(:,1),zv(:,2),zv(:,3),'o') %plots from ode113.
legend('Kepler','ode113')
hold off

%stores final position and velocity found through kepler in a row matrix
%to compare against the the final values found through ode113
x = [transpose(r) transpose(v)]; 

%zv(100,i): column 100 is the last column so it stores the final postion
%and velocity values. these values are ideally similar to the one found
%using the kepler solver in part c. i values 1-3 are postion, 4-6 velocity
for i = 1:6
    percentdiff(i) = 100.*(x(i)-zv(300,i))./x(i);
end

fprintf('a) The eccentric anomaly is: %.8f deg\n',rad2deg(E));
fprintf('b) The true anonmaly is: %.8f deg\n',rad2deg(nu))
fprintf('c) The ECI position is: [%.8f\t%.8f\t%.8f] km\n',r)
fprintf('   The ECI velocity is:  [%.8f\t%.8f\t%.8f] km/s\n',v)
fprintf('d) The differences between the postion values acquired')
fprintf(' from kepler and ode113 are: \n')
fprintf('   [%.8f\t%.8f\t%.8f] km\n',percentdiff(1),percentdiff(2),percentdiff(3))
fprintf('   The differences between the velocity values acquired')
fprintf(' from kepler and ode113 are:\n ')
fprintf('  [%.8f\t%.8f\t%.8f] km/s\n',percentdiff(4),percentdiff(5),percentdiff(6))
