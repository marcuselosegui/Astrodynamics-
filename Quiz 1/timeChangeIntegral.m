function deltat = timeChangeIntegral(f,nu1,nu2,c,N)
f = @(nu,c) ((c(1)./(1+c(2).*cos(nu))).^2)./c(3); %function created in question 1
[x,w] = GaussPointsWeights(nu1,nu2,N); %calling the gauss function givn
wt = transpose(w); %tranpsoing the gauss weights to be a row vector
fv = f(x,c); %creating a vector of the function with every x acquired from Gauss
deltat = dot(wt,fv); %multiplying the gauss weights and function to get a scalar value in seconds
end