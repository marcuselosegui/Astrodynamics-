%Question 4
clc;clear;

mu = 398600;

%orbit 1
inc = deg2rad(57);
Omega = deg2rad(60);
omega = 0;
e = 0;
nu = 0; 
Re = 6371; %radius of earth
h = 300;
r1 = Re+h;
oe1 = [r1;e;Omega;inc;omega;nu];
[rv1,vv1] = oe2rv_Elosegui_Marcus(oe1,mu);
hv1 = cross(rv1,vv1);
h1 = norm(hv1);

%orbit 2
Omega2 = 0;
omega2 = 0;
inc2 = 0;
nu2 = 0;
tauh = 23.934;
tau  = tauh*3600;
r2 = (mu*(tau/(2*pi))^2)^(1/3);
oe2 = [r2;e;Omega2;inc2;omega2;nu2];
[rv2,vv2] = oe2rv_Elosegui_Marcus(oe2,mu);
hv2 = cross(rv2,vv2);
h2 = norm(hv2);

%Transfer orbit
a = (r1+r2)/2;
e = (r2-r1)/(r2+r1);
oe3 = [a;e;Omega;inc;omega;nu];

crank = dot(hv1,hv2)/(h1*h2);

v1m = sqrt(mu/r1);
v1p = sqrt(mu)*sqrt(2/r1 - 1/a);
v2m = sqrt(mu)*sqrt(2/r2 - 1/a);
v2p = sqrt(mu/r2);

f = 0;
dv1 = sqrt(v1m^2+v1p^2-2*v1m*v1p*cos(f*crank));
dv2 = sqrt(v2m^2+v2p^2-2*v2m*v2p*cos((1-f)*crank));
dV = dv1+dv2;
   

tau_TO = 2*pi*sqrt(a^3/mu)/3600;
t = tau_TO/2;

Isp = 320;
g0 = 9.81;
m_ratio1 = exp((dv1/(g0*Isp)));
m_ratio2 = exp((dv2/(g0*Isp)));

lv = cross(hv1,hv2)/norm(cross(hv1,hv2)); %line of intersection
l = norm(lv);

u1 = cross(hv1,lv)/(h1);
u2 = -cross(hv2,lv)/(h2);
dv1v = dv1*u1;
dv2v = dv2*u2;

nu = [0:0.01:2*pi];
for i=1:length(nu)
    oe_i = [oe1(1:5); nu(i)];
    oe_ii = [oe2(1:5);nu(i)];
    oe_iii = [oe3(1:5);nu(i)];
    [RV1(i,:),VV1(i,:)] = oe2rv_Elosegui_Marcus(oe_i,mu);
    [RV2(i,:),VV2(i,:)] = oe2rv_Elosegui_Marcus(oe_ii,mu);
    [RV3(i,:),VV3(i,:)] = oe2rv_Elosegui_Marcus(oe_iii,mu);
end

plot3(RV1(:,1),RV1(:,2),RV1(:,3))
hold on
plot3(RV2(:,1),RV2(:,2),RV2(:,3))
plot3(RV3(1:end/2+1,1),RV3(1:end/2+1,2),RV3(1:end/2+1,3))
rt0 = -lv*r1;
rtf = lv*r2;
quiver3(rt0(1),rt0(2),rt0(3),dv1v(1),dv1v(2),dv1v(3),5000)
quiver3(rtf(1),rtf(2),rtf(3),dv2v(1),dv2v(2),dv2v(3),5000)

R = r2/r1;
x = [RV3(end,1),RV3(end,2),RV3(end,3)];
y = [-R*RV3(end,1),-R*RV3(end,2),-R*RV3(end,3)];
lplot = [x;y];
plot3(lplot(:,1),lplot(:,2),lplot(:,3))

fprintf('a) The magnitude of each impusle is:\n')
fprintf('   dV1 = %.8f km/s\n',dv1)
fprintf('   dV2 = %.8f km/s\n',dv2)
fprintf('b) The total dV required: %.8f km/s\n',dV)
fprintf('c) The time required: %.8f hours\n',t)
fprintf('d) The ratio mo/mf:\n')
fprintf('   for impulse 1: %.8f\n',m_ratio1)
fprintf('   for impulse 2: %.8f\n',m_ratio2)
fprintf('f) Changing the longitude of the ascending node\n')
fprintf('   does change the location on the initial orbit\n')
fprintf('   where the transfer starts\n')

