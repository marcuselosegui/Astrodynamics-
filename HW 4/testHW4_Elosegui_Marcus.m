%Test 

mu                        = 398600;
t0                        = 27*60;
tf                        = 57*60;
rv                     = [-10515.45; -5235.37; 49.17];
vv                     = [-2.10305; -4.18146; 5.563290];
[r,v,E,nu] = propagateKepler_Elosegui_Marcus(rv,vv,t0,tf,mu);
disp(r)
disp(v)
disp(E)
disp(nu)