%Question 1-8
clc; clear;

%givens
mu = 1;
r = [0; 2; 0];
v = [-1/sqrt(3); sqrt(2/3); 0];

h = cross(r,v);
e = (cross(v,h))/mu - r/norm(r);
proof = dot(h,e); 
p = (norm(h))^2 / mu;
a = p/(1-(norm(e))^2);
nu = atan2(-r(2),-r(1)) + pi;

fprintf('(a)The specific angular momentum: [%.4f\t%.4f\t%.4f]\n',h(1),h(2),h(3))
fprintf('(b)The eccentricity vector: [%.4f\t%.4f\t%.4f]\n',e(1),e(2),e(3))
fprintf('(c)h dot e = %d\n',proof)
fprintf('(d)The semi-major axis: %.4f\n',a)
fprintf('(e)The semi-latus rectum: %.4f\n',p)
fprintf('(f)The true anomaly: %.4f\n',nu)