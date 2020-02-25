function oe = rv2oe_Elosegui_Marcus(rv,vv,mu)

%-%----------------------------------------------------------------------%-
%-%------------ Code Template for Computing Orbital Elements ------------%-
%-%----------- from Cartesian Planet-Centered Inertial (PCI) ------------%-
%-%-------- and Cartesian Planet-Centered Inertial (PCI) Velocity -------%-
%-%----------------------------------------------------------------------%-
%-%--- PRIOR TO SUBMISSION THIS FUNCTION MUST BE RENAMED AS FOLLOWS: ----%-
%-%------------------------- FUNCTION HEADER LINE -----------------------%-
%-%---------- function oe = rv2oe_LastName_FirstName(rv,vv,mu) ----------%-
%-%------------------- NAME OF ACTUAL FUNCTION FILE ---------------------%-
%-%-------------------- rv2oe_LastName_FirstName.m ----------------------%-
%-%----------------------------------------------------------------------%-
%-%----------------------------------------------------------------------%-
% Inputs:                                                                %-
%    rPCI:  Cartesian planet-centered inertial (PCI) position (3 by 1)   %-
%    vPCI:  Cartesian planet-centered inertial (PCI) velocity (3 by 1)   %-
%    mu:    gravitational parameter of centrally attacting body.         %-
% Outputs:  orbital elements                                             %-
%    oe(1): semi-major axis.                                             %-
%    oe(2): eccentricity.                                                %-
%    oe(3): longitude of the ascending node (rad)                        %-
%    oe(4): inclination (rad)                                            %-
%    oe(5): argument of the periapsis (rad)                              %-
%    oe(6): true anomaly (rad)                                           %-
%-%----------------------------------------------------------------------%-
%-%----------------------------------------------------------------------%-
%-% ------------------ Format of Final Line of Code -------------------- %- 
%-% oe = [a; e; Omega; inc; omega; nu];                                  %-
%-%----------------------------------------------------------------------%-
%-%--- IMPORTANT: THE OUTPUT oe MUST BE A COLUMN VECTOR OF LENGTH SIX ---%-
r = norm(rv); %magnitude of the position vector
v = norm(vv); %magnitude of the velocity vector

hv = cross(rv,vv); %Specific angular momentum
h = norm(hv); %magnitude of the specific angular momentum
p = h^2/mu; %semi-latus rectum
ev = cross(vv,hv)./mu -rv./r; %eccentricity vector
e = norm(ev); %eccentricity 
a = p/(1-e^2); %semi-major axis

%Planet centered inertial coordinates
Ix = [1;0;0];
Iy = [0;1;0];
Iz = [0;0;1];
I = [Ix Iy Iz];

%Ascending node
nv = cross(Iz,hv);
n = norm(nv);
Omega = atan2(-dot(nv,Iy),-dot(nv,Ix)) + pi; %longitude of the ascending node

inc = atan2(dot(hv,cross(nv,Iz)),n*dot(hv,Iz)); %inclination

omega = atan2(-dot(ev,cross(hv,nv)),-h*dot(ev,nv)) + pi; %Argument of the periapsis

nu = atan2(-dot(rv,cross(hv,ev)),-h*dot(rv,ev)) + pi; %True anomaly

oe = [a; e; Omega; inc; omega; nu]; %orbital elements
end
