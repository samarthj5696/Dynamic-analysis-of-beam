% To compute time displacement graph of damped cantilever beam

% ---------------------------------------------------
%intial values
rho =1; %Density
A = 1;  %Cross sectional Area
E = 1;  %Modulos of elasticity
I = 1;  %Moment of Inertia
L = 1;  %Lenght
c1 = 0.5;c2 = 0.5;  %Damping Coefficient
Ne = 50;   %Number of elements

xx = 1/Ne:1/Ne:1;
D0 = zeros(2*Ne,1); %Intial Displacement
V0 = zeros(2*Ne,1); %Intial Velocity

%beam3 function called to find Global Mass Matrix and Stiffness Matrix
[Ma,Ka] = Beam_Matrix(rho,A,E,I,L/Ne,Ne+1,'cantilever'); 

Ca = c1*Ma+c2*Ka; %Rayleigh Damping Coefficient 
dt = 1e-4; %10^-3 Time inteval
T = 15;  %Total time
F = @(t) ExternalForce(t,Ne); %External forces 

%Newmark function called to compute the D, V, and A
[Displacement] = Newmark(Ma,Ca,Ka,F,D0,V0,dt,T); 
% -------------------------------------------------------------------------
%Plotting

plot(0:dt:T,Displacement(2*Ne-1,:),'linewidth',1)
xlabel('t')
title('Displacement of end-point')

% -------------------------------------------------------------------------


function f = ExternalForce(t,Ne)
% the external force is customizable
f = zeros(2*Ne,1);

f(2*Ne-1) =5;

end