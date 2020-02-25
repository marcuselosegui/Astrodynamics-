function [R,dVp,dVa,dV,t] = twoNImpulseOrbitTransfer(r0,rf,inc0,incf,Omega0,Omegaf,N,mu)
    
    n = 1;
    for i = 1:N
        rn(i) = r0 + (n/N)*(rf-r0);
        Omegan(i) = Omega0 + (n/N)*(Omegaf-Omega0);
        incn(i) = inc0 + (n/N)*(incf-inc0);
        n = n+1;
    end
    for ii = 1:N
        %orbit 1
        e1 = 0; %circular orbit 
        omega1 = 0; %circular orbit
        nu1 = 0; %circular orbit
        if n = 1
            r1 = r0;
            Omega1 = Omega0;
            inc1 = inc0;
        else
            r1 = rn(ii-1);
            Omega1 = Omegan(ii-1);
            inc1 = inc0(ii-1);
        end
        
        oe1 = [r1;e1;Omega1;inc1;omega1;nu1];
        [rv1,vv1] = oe2rv_Elosegui_Marcus(oe1,mu);
        hv1 = cross(rv1,vv1);
        h1 = norm(hv1);

        %orbit 2
        e1 = 0; %circular orbit 
        omega1 = 0; %circular orbit
        nu1 = 0; %circular orbit
        r2 = rn(ii);
        Omega2 = Omegan(ii);
        inc2 = incn(ii);
        oe2 = [r2;e2;Omega2;inc2;omega2;nu2];
        [rv2,vv2] = oe2rv_Elosegui_Marcus(oe2,mu);
        hv2 = cross(rv2,vv2);
        h2 = norm(hv2);
        
        %Transfer orbit
        a = (r1+r2)/2;
        e = (r2-r1)/(r2+r1);
        
        v1m = sqrt(mu/r1);
        v1p = sqrt(mu)*sqrt(2/r1 - 1/a);
        v2m = sqrt(mu)*sqrt(2/r2 - 1/a);
        v2p = sqrt(mu/r2);
        
        crank = dot(hv1,hv2)/(h1*h2);
        
        dv1 = zeros(N,1);
        dv2 = zeros(N,1);
        dv  = zeros(N,1);
        f = 0;
        dv1(ii) = sqrt(v1m^2+v1p^2-2*v1m*v1p*cos(f*crank));
        dv2(ii) = sqrt(v2m^2+v2p^2-2*v2m*v2p*cos((1-f)*crank));
        dv(ii) = dv1(ii)+dv2(ii);
        
        %for time, find tau/2 for each tranfer orbit and sum them together
    
    end
    dV = sum(dv);
end