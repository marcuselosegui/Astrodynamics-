%Task 2
clc; clear;

mu = 398600; %km^3/s
rv = [-1217.39430415697; -3091.41210822807; -6173.40732877317]; %km
vv = [9.88635815507896; -0.446121737099303; -0.890884522967222]; %km/s

oe = rv2oe_Elosegui_Marcus(rv,vv,mu); 

a = oe(1);
e = oe(2);
p = a*(1-e^2);
nu0 = oe(6);
nu_a = deg2rad(180); %apoapsis
Ni = 20; 
f = [];
t0 = timeChangeIntegral(f,nu_a,nu0,p,e,mu,Ni); %seconds
t0_m = t0/60; %minutes
tau = 2*pi*sqrt(a^3/mu); %seconds time period
tau_m = tau/60; %minutes
tn = t0+ 2*tau;
tn_m = tn/60;
n = int16((tn_m-t0_m)/5)+1; %number of increments. '+1' to account for t0

t = zeros(n,1);
t(1) = t0;
N = zeros(n,1);
N(1) = nu0;
nui = nu0;
ti = t0; %lower time bound
tii = t0 + 300; %upper time bound

for i = 2:n
   N(i) = rootFinder(f,ti,nui,tii,p,e,mu,Ni);
   t(i) = tii; %storing in time matrix
   ti = tii; %updating lower time bound
   tii = tii + 300; %updating upper time bound
   nui = N(i); %updating lower nu bound
end

%Task 3
R = zeros(n,3);
V = zeros(n,3);
for i = 1:n
    oe_i = [oe(1:5); N(i)];
    [R(i,:),V(i,:)] = oe2rv_Elosegui_Marcus(oe_i,mu);
end

%Task 4
earthSphere(50)
hold on
plot3(R(:,1),R(:,2),R(:,3),'LineWidth',2)
view(49.5,22.8)

%Task 5
lonInertial = zeros(n,1); %theta_i
lat   = zeros(n,1); %phi_i
lonEarth = zeros(n,1); %theta_e
omega_e = 2*pi/86400; %rotation rate in rad/s

x = R(:,1);
y = R(:,2);
z = R(:,3);

lonInertial = atan2(y,x);
lat = atan2(z,sqrt(x.^2+y.^2));
lonEarth = lonInertial-omega_e*t;
lonEarth = mod(lonEarth,2*pi)-pi;

earth = imread('earth.jpg');
figure(2)
clf
image('CData',earth,'XData',[-180 180],'YData',[90 -90])
hold on
plot(lonEarth*180/pi,lat*180/pi,'*');

%Task 6
%process was to find the index number n where the desired dvalues are
%wanted. lonEarth was rounded to whole number so that number like 80 or 100
%could be foun. The index n where that occurs is then used in the N matrix
%that stores the nu values. 
lonEarthR = round(rad2deg(lonEarth));
latR = round(rad2deg(lat));
n1 = find(lonEarthR == 80);
n2 = find(lonEarthR == 100);
n3 = find(lonEarthR == -100);
n4 = find(lonEarthR == -80);
deltati = timeChangeIntegral(f,N(n1),N(n2),p,e,mu,Ni);
deltatii = timeChangeIntegral(f,N(n3),N(n4),p,e,mu,Ni);

maximum = max(lat);
n_lm =  find(lat == maximum); %local maximum 1 index number
%local max 2 is found two n_lm away from n_lm since graph is symmetrical
%therefore second max = n_lm+2*n_lm = 3*n_lm
lonEarthm1 = lonEarth(n_lm)*180/pi;
latm1 = lat(n_lm)*180/pi;
lonEarthm2 = lonEarth(3*n_lm)*180/pi;
latm2 = lat(n_lm)*180/pi;

fprintf('----------------------------------------------------------\n');
fprintf('----------------------------------------------------------\n');
fprintf(' Part (a): Time spent in each segment of the orbit: [sec] \n');
fprintf(' [80 100]   deg: %16.8f                             \n',deltati);
fprintf(' [-100 -80] deg: %16.8f                             \n',deltatii);
fprintf('----------------------------------------------------------\n');
fprintf('----------------------------------------------------------\n');
fprintf(' Part (b): Both segments correspond to the apoapsis       \n');
fprintf('----------------------------------------------------------\n');
fprintf('----------------------------------------------------------\n');
fprintf(' Part (c): The geocentric latitude attains a local maximum\n');
fprintf(' at (%.8f,%.8f) and (%.8f,%.8f).        \n',lonEarthm1,latm1,lonEarthm2,latm2);
fprintf('----------------------------------------------------------\n');
fprintf('----------------------------------------------------------\n');
fprintf(' Part (d): Neither points lies within a city. The first   \n');
fprintf(' point (%.8f,%.8f) lies in Russia, in the \n',lonEarthm1,latm1);
fprintf(' middle of Siberia. The second point (%.8f,%.8f)   \n',lonEarthm2,latm2);
fprintf(' lies in the Hudson Bay in Canada.\n');
fprintf('----------------------------------------------------------\n');
fprintf('----------------------------------------------------------\n');
fprintf(' Part (e): The spacecrafts spends most of its time in the \n');
fprintf(' ranges in part (a) due to the large inclination of the orbit\n');
fprintf(' and the large eccentricity vector. The large inclination \n');
fprintf(' orients the orbit plain so that the northern hemisphere  \n');
fprintf(' views the apoapsis, which is where the spacecraft spends \n');
fprintf(' the vast majority of its time. The first figure provides \n');
fprintf(' a visulaization of this. It spends most of its time there\n'); 
fprintf(' due to the large size, i.e. large area, of the apoapsis side\n');
fprintf(' of the orbit. And larger area means larger time.\n');