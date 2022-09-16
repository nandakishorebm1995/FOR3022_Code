function [U_recon] = FML_FEM_ROM(E,K,L0,dt,tsim,Phi)
% This is a reduced-order model function for the high-fidelity model. It
% returns out the reduced-order solution

% Input:    
%           - E     : Full mass matrix
%           - K     : Full stiffness matrix
%           - L0    : Full load matrix
%           - dt    : Time step
%           - tsim  : Total simulation time
%           - Phi   : ROBs

% Output:   U_recon : Reduced-order solution
            

%% REDUCED ORDER SOLUTION

% Reduced system matrices
Er = Phi'*E*Phi; Kr = Phi'*K*Phi; Fr = Phi'*L0;
sp = size(Kr,1);
[alpha] = NewmarkIntegration3(Er,Kr,Fr,dt,sp,tsim);
U_recon = Phi*alpha;

end

