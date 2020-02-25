%Question 3
clc; clear;

R1 = 1.496e8;
R2 = 2.279e8;
mu_s = 1.327e11;

tau1 = 31.6e6; %orbital period for earth
tau2 = 59.4e6; %orbital period for mars
 
n1 = 2*pi/tau1; %mean motion of Earth
n2 = 2*pi/tau2; %mean motion of mars


t_tau = pi*sqrt((R1+R2)^3/(8*mu_s));
phase_0 = pi - n2*t_tau; 
phase_f = pi- n1*t_tau; 
phase_0p = -phase_f;
tw = -(2*phase_f)/(n2-n1);

k = 1;
if tw < 0
    if n1>n2
        tw = -(2*phase_f + 2*pi*k)/(n2-n1);
    else
        tw = -(2*phase_f - 2*pi*k)/(n2-n1);
    end
end

a = (R1+R2)/2;
dV1 = sqrt(2*mu_s/R1 - mu_s/a)-sqrt(mu_s/R1);
dV2 = sqrt(mu_s/R2) - sqrt(2*mu_s/R2 - mu_s/a);

fprintf('a) The mean motion of:\n')
fprintf('   Earth: %.8f rad/s\n',n1)
fprintf('   Mars: %.8f rad/s\n',n2)
fprintf('b) The phase angle is: %.8f deg\n',rad2deg(phase_0))
fprintf('c) The phase angle is: %.8f deg\n',rad2deg(phase_0p))
fprintf('d) The time available is: %.8f sec\n',tw)
fprintf('e) The impulses required are:\n')
fprintf('   dV1 = %.8f\n',dV1)
fprintf('   dV2 = %.8f\n',dV2)


