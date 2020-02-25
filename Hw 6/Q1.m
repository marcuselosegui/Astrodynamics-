%Question 1
clc; clear;

mu_e = 3.986e5;
mu_s = 1.327e11;
mu_m = 4.283e4;
h = 350;
Re = 6371;
r1 = Re + h;
Rm = 3389.5;
r2 = Rm + h;

R1 = 1.496e8;
R2 = 2.279e8;

vinf = sqrt(mu_s/R1)*(sqrt((2*R2)/(R1+R2))-1);
vp = sqrt(vinf^2 + 2*mu_e/r1);
vc1 = sqrt(mu_e/r1);

dVesc = vp - vc1;

e = 1 + (r1*vinf^2)/mu_e;
beta = acos(1/e);

vinf2 = sqrt(mu_s/R2)*(1-sqrt((2*R1)/(R1+R2)));
vp2 = sqrt(vinf2^2 + 2*mu_m/r2);
vc2 = sqrt(mu_m/r2);
dVcap = vp2 - vc2;

fprintf('a) The impulse required is: %.8f km/s\n',dVesc)
fprintf('b) Beta: %.8f deg\n',rad2deg(beta))
fprintf('c) The impulse required is: %.8f km/s\n',dVcap)