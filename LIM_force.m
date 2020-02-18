function [F,Xm,G,K] = LIM_force(f,p,Vs,u_0,ge,Pr,d,Wse,Kw,N1,S,I1,Ls)
%LIM_force calculates thrust-force for a given linear 
%induction motor with parameters entered as follows:
% 
% f = Supply frequency in Hertz
% p = number of poles
% Vs = 2*f*tao, synchronous velocity
% u_0 = magnetic permeability 
% ge = Air gap between primary and secondary measured in meters
% Pr = Volume-resistivity of rotor conductor outer layer
% d = Thicknes of rotor conductor outer layer
% Wse = Statorwidth
% Kw = Winding factor
% N1 = Number of windings
% S = Slipfactor
% I1 = Rated statorcurrent
% Ls = Length of stator
% 
% SYNTAX:
% F = LIM_force(f,p,Vs,u_0,ge,Pr,d,Wse,Kw,N1,S,I1,Ls)

Tao  = Ls/p; %pole-pitch
r = Kw*Tao/2;
Xm = 24*u_0*pi*f*Wse*Kw*(N1^2)*r;
G = (2*u_0*f*(Tao^2))/(pi*(Pr/d));
 

R2 = Xm/G;
K = ((1/(S*G)^2)+1)*Vs*S;
F = (3*(I1^2)*R2)/K;
end

