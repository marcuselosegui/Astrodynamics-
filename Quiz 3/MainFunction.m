%Task 1
clc;clear;

mu = 398600;
tauh = 24; %hours
taum = 24*60; %minutes
tau = tauh*3600; %seconds
e = 0.25;
inc = deg2rad(63.4); %orbital inclination 
omegap_deg = [90; 270];
omegap = deg2rad(omegap_deg); %argument of the periapsis
t0 = 0;
Omegan_deg = [0;45;90;135;180]; %long of ascending node
Omegan = deg2rad(Omegan_deg);

tinterval_min = [t0:5:taum]';
tinterval = tinterval_min*60;
tintervalh = tinterval_min/60;
n = length(tinterval);
a  = (mu*(tau/(2*pi))^2)^(1/3);
nu0 = 0; 
nomega = length(omegap);
nOmega = length(Omegan);

lonInertial = zeros(n,2*nOmega); %theta_i
lat   = zeros(n,2*nOmega); %phi_i
lonEarth = zeros(n,2*nOmega); %theta_e
omega_e = 2*pi/86400; %rotation rate in rad/s
DeltaT = zeros(2,1);
for i = 1:nomega
    for ii = 1:nOmega
        oe = [a;e;Omegan(ii);inc;omegap(i);nu0];

        [r0,v0] = oe2rv_Elosegui_Marcus(oe,mu);

        RR = zeros(n,3);
        VV = zeros(n,3);
        RR(1,1:3) = transpose(r0);
        VV(1,1:3) = transpose(v0);
        for iii = 2:n
            [RR(iii,:),VV(iii,:)] = propagateKepler_Elosegui_Marcus(transpose(RR(iii-1,1:3)),transpose(VV(iii-1,1:3)),tinterval(iii-1),tinterval(iii),mu);
        end
        
        if i == 1
            figure(1)
            earthSphere(50)
            hold on
            plot3(RR(:,1),RR(:,2),RR(:,3),'LineWidth',2)
            title(sprintf('omega = %d', omegap_deg(1)));
            
            x(:,ii) = RR(:,1);
            y(:,ii) = RR(:,2);
            z(:,ii) = RR(:,3);
            lonInertial(:,ii) = atan2(y(:,ii),x(:,ii));
            lat(:,ii) = atan2(z(:,ii),sqrt(x(:,ii).^2+y(:,ii).^2));
            lonEarth(:,ii) = lonInertial(:,ii)-omega_e*tinterval;
            lonEarth(:,ii) = mod(lonEarth(:,ii),2*pi)-pi;
            
        else
            figure(2)
            earthSphere(50)
            hold on
            plot3(RR(:,1),RR(:,2),RR(:,3),'LineWidth',2)
            title(sprintf('omega = %d', omegap_deg(2)));
            
            x(:,ii+5) = RR(:,1);
            y(:,ii+5) = RR(:,2);
            z(:,ii+5) = RR(:,3);
            lonInertial(:,ii+5) = atan2(y(:,ii+5),x(:,ii+5));
            lat(:,ii+5) = atan2(z(:,ii+5),sqrt(x(:,ii+5).^2+y(:,ii+5).^2));
            lonEarth(:,ii+5) = lonInertial(:,ii+5)-omega_e*tinterval;
            lonEarth(:,ii+5) = mod(lonEarth(:,ii+5),2*pi)-pi;
            
        end 
    end
   
    
    earth = imread('earth.jpg');
    figure(i+2)
    clf
    image('CData',earth,'XData',[-180 180],'YData',[90 -90])
    hold on
    k = 1;
    if i > 1
        i = 6;
        k = 2;
        
    end
    plot(lonEarth(:,i)*180/pi,lat(:,i)*180/pi,'*');
    plot(lonEarth(:,i+1)*180/pi,lat(:,i+1)*180/pi,'o');
    plot(lonEarth(:,i+2)*180/pi,lat(:,i+2)*180/pi,'s');
    plot(lonEarth(:,i+3)*180/pi,lat(:,i+3)*180/pi,'d');
    plot(lonEarth(:,i+4)*180/pi,lat(:,i+4)*180/pi,'v');
    title(sprintf('omega = %d', omegap_deg(k)));
    
    ff = find(lonEarth(1:end-1,i)>=lonEarth(1,i));
    hh = find(diff(ff)>1); 
    tfirst = tinterval(hh + 1);
    tsecond = tinterval(ff(end));  
    DeltaT(k) = tsecond - tfirst;
    DeltaT(k) = DeltaT(k)/3600;
end

fprintf('----------------------------------------------------------------------------\n');
fprintf('                       Task 1                                               \n');
fprintf('----------------------------------------------------------------------------\n');
fprintf('The structure of the orbit as a function of Omega changes changes the same  \n');
fprintf('way for both values of omega. The longititude over which the spacecraft orbits\n');
fprintf('changes as Omega is changed. The key difference in the structure of the orbit\n');
fprintf('is that the periapsis for omega = %.f is in the northern hemisphere while the\n',omegap_deg(1));
fprintf('periapsis for omega = %.f is in the souther hemipshere.\n',omegap_deg(2));
fprintf('----------------------------------------------------------------------------\n');
fprintf('                       Task 2                                               \n');
fprintf('----------------------------------------------------------------------------\n');
fprintf('The key difference between the orbits omega = %.f and omega = %.f is that\n  ',omegap_deg(1),omegap_deg(2));
fprintf('the figure 8''s are flipped when changing omega. The spacecraft spend the\n   ');
fprintf('majority of its time in one hemisphere for the first value of omega and in the\n');
fprintf('opposite hemisphere for the other orbit.\n');
fprintf('----------------------------------------------------------------------------\n');
fprintf('                       Task 3                                               \n');
fprintf('----------------------------------------------------------------------------\n');
fprintf('(a) for omega = %.f: \n',omegap_deg(1));
fprintf(' (i) The spacecraft spends the majority of its time in the Southern Hemisphere\n');
fprintf(' (ii) The time spent in that hemisphere is: %.8f hours\n',tauh-DeltaT(1));
fprintf('(a) for omega = %.f: \n',omegap_deg(2));
fprintf(' (i) The spacecraft spends the majority of its time in the Northern Hemisphere\n');
fprintf(' (ii) The time spent in that hemisphere is: %.8f hours\n',tauh-DeltaT(2));
fprintf(' (b) for omega = %.f, the values of Omega that place the spacecraft over the\n',omegap_deg(1));
fprintf(' following locations are:\n');
fprintf(' Australia:     220 deg\n');
fprintf(' South America: 20 deg\n');
fprintf(' (c) for omega = %.f, the values of Omega that place the spacecraft over the\n',omegap_deg(2));
fprintf(' following locations are:\n');
fprintf(' Center of Russia:        0 deg\n');
fprintf(' Center of North America: 180 deg\n');

