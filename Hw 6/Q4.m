%Question 4
clc;clear;

mu_s = 1.327e11;
tau_m = 7.6e6; %orbital period of mercury
tau = 5/2 * tau_m; %orbital period after flyby
a = ((tau/(2*pi))^2 * mu_s)^(1/3);

R1 = 1.496e8; %distance from sun to Earth
R2 = 57.91e6; %distance from sun to Mercury
R3 = 1.082e8; %Distance from sun to Venus

v_plus = sqrt(2*mu_s/R3 - mu_s/a);
Vp = sqrt(mu_s/R3);
v_minus = sqrt(2*mu_s/R3 - 2*mu_s/(R1+R3));
vinf = abs(Vp-v_minus);
delta = acos((v_plus^2-vinf^2-Vp^2)/(2*Vp*vinf));

e = 1/sin(delta/2);
mu_v = 3.248e5;
rp = ((e-1)*mu_v)/vinf^2;
Rv = 6051.8; %radius of venus
h = rp - Rv;
fprintf('a) The semi-major axis: %.8f km\n',a)
fprintf('b) The turn angle: %.8f deg\n',rad2deg(delta))
fprintf('c) The altitude of the periapsis of venus: %.8f km\n',h)
fprintf('d) The eccentricity: %.8f\n',e)
fprintf('e) yes it is possible because the hyperbolic trajectory i.e. e>1\n')