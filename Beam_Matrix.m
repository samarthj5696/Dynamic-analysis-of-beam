function [Ma,Ka] = Beam_Matrix(rho,A,E,I,Le,n,type)
% ------------------------------------------------------------------------
% Discretization by finite element method
% n is the NUMBER of NODES including the left end
% -------------------------------------------------------------------------

% Local element - Mass matrix and Stiffness matrix
Me = rho*A*Le/420*[156,    22*Le,   54,     -13*Le;
                   22*Le,  4*Le^2,  13*Le,  -3*Le^2;
                   54,     13*Le,   156,    -22*Le;
                   -13*Le, -3*Le^2, -22*Le, 4*Le^2];
                   
Ke = E*I/Le^3*[12,    6*Le,    -12,    6*Le;
               6*Le,  4*Le^2,  -6*Le,  2*Le^2;
               -12,   -6*Le,   12,     -6*Le;
               6*Le,  2*Le^2,  -6*Le,  4*Le^2];
% -------------------------------------------------------------------------
% Global stiffness matrix
Ma = zeros(2*n,2*n);
Ka = zeros(2*n,2*n);
for i = 1:2:2*n-3
    Ma(i:i+3,i:i+3) = Ma(i:i+3,i:i+3) + Me;
    Ka(i:i+3,i:i+3) = Ka(i:i+3,i:i+3) + Ke;
end 
% -------------------------------------------------------------------------
% Boundary Conditions
switch type
    case 'cantilever'
        % the left end is clamped !
        Ma(:,1:2) = [];Ma(1:2,:) = [];
        Ka(:,1:2) = [];Ka(1:2,:) = [];
    case 'simply-supported'
        % simply supported at two ends
        Ma(:,[1,end-1]) = [];Ma([1,end-1],:) = [];
        Ka(:,[1,end-1]) = [];Ka([1,end-1],:) = [];
   
end
%--------------------------------------------------------------------------