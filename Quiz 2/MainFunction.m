clc;clear;

%givens
mu = 398600;
hp = 598.77493602896; %periapsis altitude
ha = 39952.9133996976; %apoapsis altitude
re = 6378.145; %radius of earth
nu1deg = 90; 
nu2deg = 270;

rp = hp+re; %periapsis radius
ra = ha+re; %apoapsis radius
nu1 = deg2rad(nu1deg);
nu2 = deg2rad(nu2deg);

e = (ra-rp)/(ra+rp); %eccentricity 
a = (rp+ra)/2; %semi-major axis
p = a*(1-e^2); %semi-latus rectum
h = sqrt(p*mu); %specific angular momentum
c = [p;e;h]; %vector c

tau = 2*pi*sqrt(a^3/mu); %orbital period
tauhr = tau/3600; %orbital period converted to hours

N = [10; 15; 20; 25];
f=[];
deltat = [];
for i = 1:length(N)
deltati = timeChangeIntegral(f,nu1,nu2,p,e,mu,N(i));
deltat = [deltat; deltati]; %time interval stored in a column matrix
end
deltat1 = deltat/3600; %time interval converted into hours

fprintf('----------------------------------------------------------\n');
fprintf('----------------------------------------------------------\n');
fprintf(' Part (a): Quantities for Orbital Period & Function \n');
fprintf('----------------------------------------------------------\n');
fprintf('----------------------------------------------------------\n');
fprintf('semi-latus rectum (c1): [km]:\t\t %16.8f\n',p);
fprintf('----------------------------------------------------------\n');
fprintf('Eccentricity (c2): [no units]:\t\t %16.8f\n',e);
fprintf('----------------------------------------------------------\n');
fprintf('specific angular momentum (c3): [km/s]:\t\t %16.8f\n',h);
fprintf('----------------------------------------------------------\n');
fprintf('----------------------------------------------------------\n');
fprintf(' Part (b): Orbital Period \n');
fprintf('----------------------------------------------------------\n');

fprintf('----------------------------------------------------------\n');
fprintf('Orbital period [hours]:\t\t\t\t\t %16.8f\n',tauhr);
fprintf('----------------------------------------------------------\n');
fprintf('----------------------------------------------------------\n');
fprintf(' Part (c): Time change from nu1=%i deg to nu2=%i deg \n',nu1deg,nu2deg)
;
fprintf('----------------------------------------------------------\n');
fprintf('----------------------------------------------------------\n');
for i=1:length(N)
    fprintf('Time elapsed (%i,%i) deg [hours] (N=%i):%16.8f\n',nu1deg,nu2deg,N(i),deltat1(i));
end
fprintf('----------------------------------------------------------\n');
fprintf('----------------------------------------------------------\n');
fprintf(' Part (d): Time change from nu2=%i deg to nu1=%i deg \n',nu2deg,nu1deg);
fprintf('----------------------------------------------------------\n');
fprintf('----------------------------------------------------------\n');
for i=1:length(N)
fprintf('Time elapsed (%i,%i) deg [hours] (N=%i):%15.8f\n',nu2deg,nu2deg+180,N(i),deltat2(i));
end
fprintf('----------------------------------------------------------\n');
fprintf('----------------------------------------------------------\n');
fprintf(' Part (e): An interesting thing to note about the orbital \n');
fprintf(' period calculated in part (b) is that it is not much larger\n');
fprintf(' than the time time intervals calculated in part (c) which \n');
fprintf(' went only from 90 to 270 degrees, or half of the orbit. \n');
fprintf(' This shows just how much larger the apoapis radius is  \n');
fprintf(' than the periapsis radius. Since equal area equals equal \n');
fprintf('time, this shows that a very large percentage of a spacraft \n');
fprintf(' in orbit would spend it in the apoapsis area of the orbit\n');
fprintf('----------------------------------------------------------\n');
fprintf('----------------------------------------------------------\n');