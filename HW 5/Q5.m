%Question 5
clc;clear;

mu = 398600;
Re = 6371;

%orbit 1
h = 350;
r1 = Re+h;
inc1 = deg2rad(28);
e1 = 0;
omega1 = 0;
Omega1 = 0;
nu1 = 0;
oe1 = [r1;e1;Omega1;inc1;omega1;nu1];
[rv1,vv1] = oe2rv_Elosegui_Marcus(oe1,mu);
hv1 = cross(rv1,vv1);
h1 = norm(hv1);

%orbit 2
r2 = 26558;
inc2 = deg2rad(55);
e2 = 0;
omega2 = 0;
Omega2 = Omega1;
nu2 = 0;
oe2 = [r2;e2;Omega2;inc2;omega2;nu2];
[rv2,vv2] = oe2rv_Elosegui_Marcus(oe2,mu);
hv2 = cross(rv2,vv2);
h2 = norm(hv2);

lv = cross(hv1,hv2)/norm(cross(hv1,hv2)); %line of intersection

%Transfer orbit
a = (r1+r2)/2;
e = (r2-r1)/(r2+r1);
nu_a = pi; %nu at apoapsis where 2nd impulse occurs
oe3 = [a;e;Omega1;inc1;omega1;nu_a];
[rv3,vv3] = oe2rv_Elosegui_Marcus(oe3,mu);

v1m = sqrt(mu/r1);
v1p = sqrt(mu)*sqrt(2/r1 - 1/a);
v2m = sqrt(mu)*sqrt(2/r2 - 1/a);
v2p = sqrt(mu/r2);

crank = dot(hv1,hv2)/(h1*h2);

f = 0;
dv1 = sqrt(v1m^2+v1p^2-2*v1m*v1p*cos(f*crank));
dv2 = sqrt(v2m^2+v2p^2-2*v2m*v2p*cos((1-f)*crank));
dV = dv1+dv2;

u1 = cross(hv1,lv)/(h1);
u2 = -cross(hv2,lv)/(h2);
dv1v = dv1*u1;
dv2v = dv2*u2;


tau_TO = 2*pi*sqrt(a^3/mu)/3600;
t = tau_TO/2;


nu = [0:0.01:2*pi];
for i=1:length(nu)
    oe_i = [oe1(1:5); nu(i)];
    oe_ii = [oe2(1:5);nu(i)];
    oe_iii = [oe3(1:5);nu(i)];
    [RV1(i,:),VV1(i,:)] = oe2rv_Elosegui_Marcus(oe_i,mu);
    [RV2(i,:),VV2(i,:)] = oe2rv_Elosegui_Marcus(oe_ii,mu);
    [RV3(i,:),VV3(i,:)] = oe2rv_Elosegui_Marcus(oe_iii,mu);
end

R = r2/r1;
plot3(RV1(:,1),RV1(:,2),RV1(:,3))
hold on   
plot3(RV2(:,1),RV2(:,2),RV2(:,3))
plot3(RV3(1:end/2+1,1),RV3(1:end/2+1,2),RV3(1:end/2+1,3))

rt0 = lv*r1;
rtf = lv*-r2;
quiver3(rt0(1),rt0(2),rt0(3),dv1v(1),dv1v(2),dv1v(3),5000)
quiver3(rtf(1),rtf(2),rtf(3),dv2v(1),dv2v(2),dv2v(3),5000)


x = [RV3(end,1),RV3(end,2),RV3(end,3)];
y = [-R*RV3(end,1),-R*RV3(end,2),-R*RV3(end,3)];
lplot = [x;y];
plot3(lplot(:,1),lplot(:,2),lplot(:,3))

fprintf('a) The line of intersection: [%d %d %d]\n',lv)
fprintf('b) The postions of the spacecraft where the impulse is applied:\n')
fprintf('   r1: [%d %d %d] km\n',rv1)
fprintf('   r2: [%d %d %d] km\n',rv2)
fprintf('c) The total dV required: %.8f km/s\n',dV)
fprintf('d) The time required: %.8f hours\n',t)
fprintf('e) The eccentricity of the TO: %.8f\n',e)
fprintf('d) The cranking angle: %.8f deg\n',rad2deg(crank))
