function deltat = timeChangeIntegral(f,nu1,nu2,p,e,mu,N)
f = @(nu,p,e,mu) sqrt(p.^3./mu)./(1.+e.*cos(nu)).^2; %function created in question 1
[x,w] = GaussPointsWeights(nu1,nu2,N); %calling the gauss function givn
wt = transpose(w); %tranpsoing the gauss weights to be a row vector
fv = f(x,p,e,mu); %creating a vector of the function with every x acquired from Gauss
deltat = dot(wt,fv); %multiplying the gauss weights and function to get a scalar value in seconds
end