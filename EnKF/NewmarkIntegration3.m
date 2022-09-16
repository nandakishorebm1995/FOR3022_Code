function [ukp1,udkp1,uddkp1] = NewmarkIntegration3(M,K,fr,dt,sp,k,uk,udk,uddk)

% This function implements the implicit Newmark-beta time integration
% method. It is an explicit integration scheme. Here, alpha and beta are
% the constants that are set to 0.25 and 0.5 respectively, such that the
% scheme remains unconditionally stable. 
% Input: 
%             M    - Mass matrix
%             K    - Stiffness matrix
%             fr   - Load matrix
%             dt   - Time step 
%             sp   - No. of rows of the output matrix
%             tsim - Total simulation time
%
% Output:
%             dircoeff - Output matrix
%%
beta = 0.25;  % Newmark parameter: 0<= beta <=0.25
gamma = 0.5;  % Newmark numerical damping parameter:  gamma = 0.5
 
ukp1 = uk;
udkp1 = udk;
uddkp1 = uddk;

uk   = ukp1;
udk  = udkp1;
uddk = uddkp1;

fr(:,size(fr,2)+1)=zeros(sp,1);

% STEP 1: Calculation of the predictors
utilde_kp1 = uk + udk*dt + (uddk*(0.5 - beta)*dt^2);
utilde_d_kp1 = udk + uddk*(1-gamma)*dt;

% STEP 2: Solution of the linear problem
T1 = (M + (K*beta*dt^2));
T2 = fr(:,k) - (K*utilde_kp1);
uddkp1 = T1\T2; 

% STEP 3: Calculation of correctors
udkp1 = utilde_d_kp1 + (uddkp1*gamma*dt);
ukp1 = utilde_kp1 + (uddkp1*beta*dt^2);

end

