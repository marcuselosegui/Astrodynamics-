function [nu,E] = KeplerSolver_Elosegui_Marcus(nu0,t0,t,a,e,mu)


E0 = 2*atan2(sqrt(1-e)*sin(nu0/2),sqrt(1+e)*cos(nu0/2));
tp = t0 - sqrt(a^3/mu)*(E0-e*sin(E0));
tau = 2*pi*sqrt(a^3/mu);
k = floor((t-tp)/tau);
C = sqrt(mu/a^3)*(t-t0)- 2*pi*k+(E0-e*sin(E0));

M0 = E0-e*sin(E0);

E_ki = M0;
for i = 1:100
    E_kii = e*sin(E_ki)+C;
    E_ki = E_kii;
end

E = E_kii;

nu = 2*atan2(sqrt(1+e)*sin(E/2),sqrt(1-e)*cos(E/2));

end