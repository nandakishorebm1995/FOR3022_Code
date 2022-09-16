clear;
clc;
close all;
format long g;

%% SINEBURST EXCITATION

fex = 120000; T=1/fex; n=5; dt=T/20; tsim=5*n*T;
Vamp = 10e-9;
UU = zeros(1,500);
for i=1:100
    t=i*dt;
    UU(i)=Vamp*sin(2*pi*fex*t)*sin(pi*fex/n*t)^2;
end

%% PARAMETER SPACE AND DEFINITION

de = 50e6:5e8:50e8;         % Stiffness of the damage
dp = 0.02:0.01:0.11;        % X-position of the damage
ds = 0.002:0.001:0.011;     % Length of the damage

% Limits of the parameters
p1 = [50e6; 50e8];
p2 = [0.02; 0.11];
p3 = [0.002; 0.011];

[DE,DP,DS] = ndgrid(de,dp,ds);
paramstable = table(DE(:),DP(:),DS(:));
params = table2array(paramstable);
clear paramstable;

%% PARAMETER INITIALIZATION

ridx = randi(1000,1,1);
ridx2 = randi(16,1,1);
 
initParam = [DE(ridx),DP(ridx),DS(ridx)];

Damage_E = initParam(1);
Damage_pos = [initParam(2),1.24e-03];
Damage_size = [initParam(3),1.2e-04];

Ne = 6;
rid = randi(1000,Ne,1);

% PARAMETER ENSEMBLE INITIATION
Parameter_ensemble_0 = [25e8, 0.07, 0.01;
                        05e8, 0.06, 0.006;
                        35e8, 0.05, 0.005;
                        10e8, 0.06, 0.008;
                        18e8, 0.04, 0.004;
                        28e8, 0.07, 0.007];

%% TEST PARAMETERS
% Measurement 1: [40e8, 0.07, 0.005];
% Measurement 2: [35e8, 0.04, 0.007];
% Measurement 3: [30e8, 0.05, 0.005];
% Measurement 4: [14e8, 0.03, 0.003];
% Measurement 5: [08e8, 0.07, 0.01];
% Measurement 6: [10e8, 0.05, 0.005];
% Measurement 7: [20e8, 0.05, 0.008];

% MEASUREMENT DATA WITH SCALING 
MD = load('MD5.mat'); 
MD = abs(MD.ROS);
Observation_ensemble_q = MD(15648,1:500)'.*ones(500,Ne);
Observation_ensemble_q_scaled = rescale(Observation_ensemble_q);
samples = [0,0,0];

%% LOADING GLOBAL REDUCED-ORDER BASES (ROBs)

warn=warning('query','all');
id=warn.identifier;
warning('off',id);

% Generated ROBs
Phi = load('AdaptivePOD_550.mat'); Phi = Phi.V;
Observation_ensemble = zeros(500,Ne);
Solution_ensemble = zeros(79266,Ne);
tic; 

%% ENSEMBLE KALMAN FILTER PROCEDURE WITH FEATURE SCALING

for k = 190:2:500 % FREQUENCY OF UPDATE
    fprintf ('k = %d \n',k);
    
    for j = 1:Ne
        fprintf ('j = %d \n',j);
        [E,K,L0] = FML_FEM_MATRIX_EXTRACT(Parameter_ensemble_0(j,1),[Parameter_ensemble_0(j,2),1.24e-03],[Parameter_ensemble_0(j,3),1.2e-04],UU);
        ROS = FML_FEM_ROM(E,K,L0,dt,tsim,Phi);
        Observation_ensemble(:,j) = abs(ROS(15648,1:500))'; 
        Solution_ensemble(:,j) = abs(ROS(:,k)); % 
    end
    
    % FEATURE SCALING
    Observation_ensemble_scaled = rescale(Observation_ensemble);
    Parameter_ensemble_0_scaled(:,1) = rescale([Parameter_ensemble_0(:,1);p1]);
    Parameter_ensemble_0_scaled(:,2) = rescale([Parameter_ensemble_0(:,2);p2]);
    Parameter_ensemble_0_scaled(:,3) = rescale([Parameter_ensemble_0(:,3);p3]);
    Parameter_ensemble_0_scaled = Parameter_ensemble_0_scaled(1:Ne,:);
    Solution_ensemble_scaled = rescale(Solution_ensemble);
    variance = rescale([3e-20;Observation_ensemble(:,1)]);
    
    % PREDICTION STAGE
    mean_Observation_ensemble = mean(Observation_ensemble,2);
    mean_Observation_ensemble_to_be_scaled = rescale([Observation_ensemble, mean_Observation_ensemble]);
    mean_Observation_ensemble_scaled = mean_Observation_ensemble_to_be_scaled(:,end);
    
    mean_Parameter_ensemble_0 = mean(Parameter_ensemble_0', 2);
    mean_Parameter_ensemble_0_to_be_scaled(1,:) = rescale([Parameter_ensemble_0(:,1)', p1(1), p1(2), mean_Parameter_ensemble_0(1)]);
    mean_Parameter_ensemble_0_to_be_scaled(2,:) = rescale([Parameter_ensemble_0(:,2)', p2(1), p2(2), mean_Parameter_ensemble_0(2)]);
    mean_Parameter_ensemble_0_to_be_scaled(3,:) = rescale([Parameter_ensemble_0(:,3)', p3(1), p3(2), mean_Parameter_ensemble_0(3)]);
    mean_Parameter_ensemble_0_scaled = mean_Parameter_ensemble_0_to_be_scaled(:,end);
    
    mean_Solution_ensemble = mean(Solution_ensemble,2);
    mean_Solution_ensemble_to_be_scaled = rescale([Solution_ensemble, mean_Solution_ensemble]);
    mean_Solution_ensemble_scaled = mean_Solution_ensemble_to_be_scaled(:,end);
    
    % COMPUTATION OF COVARIANCE MATRICES
    omo = Observation_ensemble_scaled(k,:) - ones(1,Ne)*mean_Observation_ensemble_scaled(k);    
    Cshsh = (omo*omo')/(Ne-1);    
    Cmuhsh = ((Parameter_ensemble_0_scaled' - ones(3,Ne).*mean_Parameter_ensemble_0_scaled)*omo')/(Ne-1);
    Cuhsh = ((Solution_ensemble_scaled - (mean_Solution_ensemble_scaled.*ones(79266,Ne)))*omo')/(Ne-1);
    
    % ANALYSIS STAGE 
    Gamma = variance(1);   
    denom = Gamma + Cshsh;   
    sqmsh = Observation_ensemble_q_scaled(k,:) - Observation_ensemble_scaled(k,:);   
    right = sqmsh./denom;   
    rright = [Cmuhsh; Cuhsh]*right;    
    Parameter_ensemble_0_new = rright(1:3,:)';
   
    Parameter_ensemble_0_new = Parameter_ensemble_0_scaled + Parameter_ensemble_0_new;
    Parameter_ensemble_0_rescaled(:,1) = round(rescale([Parameter_ensemble_0_new(:,1); 0; 1],p1(1),p1(2)));
    Parameter_ensemble_0_rescaled(:,2) = round(rescale([Parameter_ensemble_0_new(:,2); 0; 1],p2(1),p2(2)),3);
    Parameter_ensemble_0_rescaled(:,3) = round(rescale([Parameter_ensemble_0_new(:,3); 0; 1],p3(1),p3(2)),3);    
    Parameter_ensemble_0 = Parameter_ensemble_0_rescaled(1:Ne,:);
    
    % COLLECT THE SAMPLES
    samples = [samples; Parameter_ensemble_0];
    
    vars = {'Parameter_ensemble_0_scaled', 'Parameter_ensemble_0_rescaled'};
    clear(vars{:}); 
    
 end

toc;





