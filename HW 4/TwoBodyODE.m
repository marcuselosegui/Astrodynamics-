function drdt = TwoBodyODE(t,z)
mu = 398600;
r = sqrt(z(1)^2 + z(2)^2 + z(3)^2);
G = -(mu/r^3);
drdt = [z(4);z(5);z(6);G*z(1);G*z(2);G*z(3)];
end