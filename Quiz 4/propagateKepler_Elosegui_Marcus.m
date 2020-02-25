function [R,V,E,nu] = propagateKepler_Elosegui_Marcus(rv,vv,t0,t,mu)

oe = rv2oe_Elosegui_Marcus(rv,vv,mu);
a = oe(1);
e = oe(2);
nu0 = oe(6);


[nu,E] = KeplerSolver_Elosegui_Marcus(nu0,t0,t,a,e,mu);
% E = E;

oe_f = [oe(1:5); nu];
[R,V] = oe2rv_Elosegui_Marcus(oe_f,mu);
end
