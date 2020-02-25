function nu  = rootFinder(f,t0,nu0,t,p,e,mu,N)
tInterval = t-t0;

%Newtons Method
nu_k = nu0 + deg2rad(5);

for i = 1:N
    g_nu = timeChangeIntegral(f,nu0,nu_k,p,e,mu,N) - tInterval; %given function
    dg_nu = sqrt(p^3/mu)/((1+e*cos(nu_k))^2); %derivative of given function
    nu_kk = nu_k - (g_nu / dg_nu); %Newtons method
    nu_k = nu_kk;
end

nu = nu_k;

