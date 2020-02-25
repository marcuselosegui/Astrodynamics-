function F = root2d(R)
S = [2;5;10;11;12;15];
F(1) =  (sqrt(2*R/(1+R))-1) + sqrt(1/R).*(1-sqrt(2/(1+R))); %Hohmann
F(2) = (sqrt(2)-1)*(1+sqrt(1/R)); %bi=parabolic transfer
end