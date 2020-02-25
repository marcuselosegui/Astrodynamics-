function [rv,vv]  = oe2rv_Elosegui_Marcus(oe,mu)

%-%----------------------------------------------------------------------%-
%-%-------------------- Code Template for Computing ---------------------%-
%-%--------- Cartesian Planet-Centered Inertial (PCI) Position ----------%-
%-%------- and Cartesian Planet-Centered Inertial (PCI) Velocity --------%-
%-%------------------------ from Orbital Elements -----------------------%-
%-%----------------------------------------------------------------------%-
%-%--- PRIOR TO SUBMISSION THIS FUNCTION MUST BE RENAMED AS FOLLOWS: ----%-
%-%------------------------- FUNCTION HEADER LINE -----------------------%-
%-%-------- function [rPCI,vPCI] = oe2rv_LastName_FirstName(oe,mu) ------%-
%-%------------------- NAME OF ACTUAL FUNCTION FILE ---------------------%-
%-%-------------------- oe2rv_LastName_FirstName.m ----------------------%-
%-%----------------------------------------------------------------------%-
%-%----------------------------------------------------------------------%-
%-% Input:  orbital elements             (6 by 1 column vector)          %-
%-%   oe(1): Semi-major axis.                                            %-
%-%   oe(2): Eccentricity.                                               %-
%-%   oe(3): Longitude of the ascending node (rad)                       %-
%-%   oe(4): Inclination (rad)                                           %-
%-%   oe(5): Argument of the periapsis (rad)                             %-
%-%   oe(6): True anomaly (rad)                                          %-
%-%   mu:    Planet gravitational parameter     (scalar)                 %-
%-% Outputs:                                                             %-
%-%   rPCI:  Planet-Centered Inertial (PCI) Cartesian position           %-
%-%          (3 by 1 column vector)                                      %-
%-%   vPCI:  Planet-Centered Inertial (PCI) Cartesian inertial velocity  %-
%-%          (3 by 1 column vector)                                      %-
%-%----------------------------------------------------------------------%-

%-%----------------------------------------------------------------------%-
%-% ---------- Final Two Lines of Code (UNCOMMENT)                       %- 
%-% rPCI = < PUT YOUR FINAL PCI POSITION HERE >                          %-
%-% vPCI = < PUT YOUR FINAL PCI INERTIAL VELOCITY HERE >                 %-
%-%----------------------------------------------------------------------%-

%Defining orbital elements
a = oe(1);
e = oe(2);
Omega = oe(3);
inc = oe(4);
omega = oe(5);
nu = oe(6);

p = a*(1-e^2); %Semi-latus rectum
r = p/(1+e*cos(nu)); %position
rv_P = [r*cos(nu); r*sin(nu); 0]; %postion in P

vv_P = sqrt(mu/p)*[-sin(nu); e+cos(nu); 0]; %velocity in P

%Transformation matrices
T_N_I = [cos(Omega) -sin(Omega) 0; sin(Omega) cos(Omega) 0; 0 0 1];
T_Q_N = [1 0 0; 0 cos(inc) -sin(inc); 0 sin(inc) cos(inc)];
T_P_Q = [cos(omega) -sin(omega) 0; sin(omega) cos(omega) 0; 0 0 1];
T_P_I = T_N_I*T_Q_N*T_P_Q;

rv = T_P_I*rv_P; %Position in PCI
vv = T_P_I*vv_P; %Velocity in PCI
end
