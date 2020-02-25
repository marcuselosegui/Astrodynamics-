%Question 2
clc;clear;



fun = @(R) (sqrt(2*R/(1+R))-1) + sqrt(1/R)*(1-sqrt(2/(1+R))) - (sqrt(2)-1)*(1+sqrt(1/R));
R0 = 1;
Ra = fsolve(fun,R0);

S = [2;5;10;11;12;15];
for i = 1:length(S)
    funb = @(R) (sqrt(2*R/(1+R))-1) + sqrt(1/R)*(1-sqrt(2/(1+R))) - ((sqrt(2*R*S(i)/(1+R*S(i)))-1) + sqrt(1/(R*S(i)))*(sqrt(2/(1+S(i)))-sqrt(2/(1+R*S(i)))) +(sqrt(2*S(i)/(R+R*S(i)))-sqrt(1/R)));
    Rb(i) = fsolve(funb,R0);
end

plot(S,Rb)

fprintf('a) The value of R where the total impulseHohmann transfer is\n')
fprintf('   the same as the bi-parabolic transfer: %.8f\n', Ra)
fprintf('b) The value of R where the total impulseHohmann transfer is\n')
fprintf('   the same as the bi-elliptic transfer for:\n')
for ii = 1:length(S)
    fprintf('   S = %.d: %.8f\n',S(ii),Rb(ii))
end
