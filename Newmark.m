function [Displacement] = Newmark(M,C,K,F,D0,V0,dt,T)
                                          
% ------------------------------------------------------------------------
% Newmark method with Average acce leration method
    Beta = 1/4;
    Gamma = 1/2;

% integration constant
c1 = 1/Beta/dt^2;
c2 = Gamma/Beta/dt;
c3 = 1/Beta/dt;
c4 = 1/2/Beta-1;
c5 = Gamma/Beta-1;
c6 = (Gamma/2/Beta-1)*dt;
c7 = (1-Gamma)*dt;
c8 = Gamma*dt;
% -------------------------------------------------------------------------
A0 = M\(F(0)-K*D0-C*V0);     % initial acceleration
n = T/dt+1; %time increament
m = length(D0);
Displacement = zeros(m,n);
V = zeros(m,n);
A = zeros(m,n);
Displacement(:,1) = D0;
V(:,1) = V0;
A(:,1) = A0;
% -------------------------------------------------------------------------
Kbar = c1*M+c2*C+K;       
for i = 1:n-1
    Da = Displacement(:,i); %displacement at time i
    Va = V(:,i);    %velocity at time i
    Aa = A(:,i);    %acceleration at time i
    Fbar = F(i*dt)+M*(c1*Da+c3*Va+c4*Aa)+C*(c2*Da+c5*Va+c6*Aa);
    Displacement(:,i+1) = Kbar\Fbar;
    A(:,i+1) = c1*(Displacement(:,i+1)-Da)-c3*Va-c4*Aa;
    V(:,i+1) = Va+c7*Aa+c8*A(:,i+1);
end
end