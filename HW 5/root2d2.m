function F = root2d2(R)
S = [2;5;10;11;12;15];
F(1) =  (sqrt(2*R/(1+R))-1) + sqrt(1/R).*(1-sqrt(2/(1+R))); %Hohmann
F(2) = (sqrt(2.*R.*S./(1+R.*S))-1) + sqrt(1./(R.*S)).*(sqrt(2./(1+S))-sqrt(2./(1+R.*S))) +(sqrt(2*S./(R+R.*S))-sqrt(1./R));
end