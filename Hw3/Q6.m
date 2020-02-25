%Question 6
clc; clear;
%Givens
oe = [1.2774; 0.2512; 5.4978; 0.3155; 1.8654; 5.9156];
mu = 1;

[rv,vv]  = oe2rv_Elosegui_Marcus(oe,mu);

fprintf('The position of the spacecraft: [%.4f\t%.4f\t%.4f] AU\n',rv);
fprintf('The velocity of the spacecraft: [%.4f\t%.4f\t%.4f] AU/TU\n',vv);
    
    
    
    