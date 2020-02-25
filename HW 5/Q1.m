%Question 1
clc; clear;
R = (1:20);
R2 = (10:16)
S = [2;5;10;11;12;15];
%all of these are normalized by Vc1
dVH =  (sqrt(2.*R./(1+R))-1) + sqrt(1./R).*(1-sqrt(2./(1+R))); %Hohmann
%Bi-elliptic transfer
dVBE = (sqrt(2.*R.*S./(1+R.*S))-1) + sqrt(1./(R.*S)).*(sqrt(2./(1+S))-sqrt(2./(1+R.*S))) +(sqrt(2*S./(R+R.*S))-sqrt(1./R));
%bi=parabolic transfer
dVBP = (sqrt(2)-1).*(1+sqrt(1./R));

figure(1)
plot(R,dVH,R,dVBE,R,dVBP)
figure(2)
plot(R,dVH,R,dVBE,R,dVBP)
xlim([10 16])
ylim([0.51 0.55])