%Question 3
clc;clear;

mu = 1;
vc1 = 1;
v1m = 1; %v1 minus
r1  = 1;
v2p = 0.5; %v1 plus
r2  = 4;
R = r2/r1;

e1 = 0;
a1 = r1;
inc = 0;
Omega = 0;
omega = 0;
nu = [0:0.01:2*pi];
a2 = r2;
e2 = 0;

dVH1 = ((sqrt(2.*R./(1+R))-1)); %Hohmann
dVH2 = (sqrt(1./R).*(1-sqrt(2./(1+R))));
a = (r1+r2)/2;
e = (r2-r1)/(r1+r2);
oe = [a1;e1;Omega;inc;omega;nu(1)];
oe2 = [a2;e2;Omega;inc;omega;nu(1)];
oe3 = [a;e;Omega;inc;omega;nu(1)];
for i=1:length(nu)
    oe_i = [oe(1:5); nu(i)];
    oe_ii = [oe2(1:5);nu(i)];
    oe_iii = [oe3(1:5);nu(i)];
    [rv(i,:),vv(i,:)] = oe2rv_Elosegui_Marcus(oe_i,mu);
    [rv2(i,:),vv2(i,:)] = oe2rv_Elosegui_Marcus(oe_ii,mu);
    [rv3(i,:),vv3(i,:)] = oe2rv_Elosegui_Marcus(oe_iii,mu);
end

n = length(rv3(1:end/2,1));
plot3(rv(:,1),rv(:,2),rv(:,3))
hold on
plot3(rv2(:,1),rv2(:,2),rv2(:,3))
plot3(rv3(1:end/2,1),rv3(1:end/2,2),rv3(1:end/2,3))
quiver(rv(1,1),rv(1,2),rv(1,3),dVH1)
quiver(rv3(n,1),rv3(n,2),rv3(n,3),dVH2)


