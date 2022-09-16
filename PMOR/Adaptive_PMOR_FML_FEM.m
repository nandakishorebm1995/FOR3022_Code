clear;
clc;
close all;

%% SINEBURST EXCITATION

fex = 120000; T=1/fex; n=5; dt=T/20; tsim=5*n*T;
Vamp = 10e-9;
UU = zeros(1,500);
for i=1:100
    t=i*dt;
    UU(i)=Vamp*sin(2*pi*fex*t)*sin(pi*fex/n*t)^2;
end

%% PARAMETER SPACE AND DEFINITION

de = 50e6:5e8:50e8;         % Stiffness of the damage#+
dp = 0.02:0.01:0.11;        % X-position of the damage
ds = 0.002:0.001:0.011;     % Size of the damage

[DE,DP,DS] = ndgrid(de,dp,ds);
paramstable = table(DE(:),DP(:),DS(:));
params = table2array(paramstable);
clear paramstable;

%%
ridx = randi(1000,1,1);
ridx2 = randi(17,1,1);
DZ = [0, 0.00012, 0.000245, 0.00037, 0.000495, 0.000620, 0.00074, ...
      0.000865, 0.00099, 0.001115, 0.00124, 0.00136, 0.001485, 0.00161, ...
      0.001735, 0.00186, 0.00198]; % Total 17 z-coords for 16 layers 
initParam = [DE(ridx),DP(ridx),DS(ridx), DZ(ridx2)];

% initParam(4) represents the position of the damage on the z-axis
% DH denotes the height of the damage. Since the steel and CFRP lamina has
% different heights, the following if loop is used.

if initParam(4)==0 || initParam(4)==0.00062 || initParam(4)==0.00124 || initParam(4)==0.00186
    DH = 1.2e-04;
else 
    DH = 1.25e-04;
end

Damage_E = initParam(1);
Damage_pos = [initParam(2),initParam(4)];
Damage_size = [initParam(3),DH];

%% COMPUTATION OF HIGH DIMENSIONAL MODEL (HDM)

warn=warning('query','all');
id=warn.identifier;
warning('off',id);

% EXTRACTING SYSTEM MATRICES FOR THE SPECIFIED PARAMETER 
[E,K,L0] = FML_FEM_MATRIX_EXTRACT(Damage_E,Damage_pos,Damage_size,UU);

% GENERATED GLOBAL REDUCED-ORDER BASES (ROBs)
V = load('AdaptivePOD_10.mat'); V = V.V;

% HIGH-FIDELITY SOLUTION
tic;HDM = FML_FEM_HDM(Damage_E,Damage_pos,Damage_size);toc;

% REDUCED-ORDER SOLUTION
tic;ROS = FML_FEM_ROM(E,K,L0,dt,tsim,V);toc;

% PLOTTING HIFI AND REDUCED-ORDER SOLUTIONS
figure; plot(linspace(0,tsim,496),HDM(15648,5:end),'linewidth',2);
hold on; plot(linspace(0,tsim,499),ROS(15648,2:end),'linewidth',2);

%% PMOR PARAMETER SPECIFICATION AND SVD TO OBTAIN THE FIRST MODE
[u,s,v] = svd(HDM,'econ');

% MAXIMUM NUMBER OF GLOBAL MODES TO BE EXTRACTED
Nitermax = 500;

N_c = 20;

V = u(:,1);

%% ADAPTIVE PMOR GREEDY PROCEDURE ACCORDING TO ALGORITHM 2 IN 
% "Parametric Model Order Reduction of Guided Ultrasonic Wave Propagation 
% in Fiber Metal Laminates with Damage"

muNiters = zeros(Nitermax,3);
apostError = zeros(1,Nitermax);
pointError = zeros(1,Nitermax);

for i = 2:Nitermax
    
    disp(i)
    fprintf('------------\n');
    Ncinit = 10;
    Nc0 = params(randperm(size(params,1),Ncinit),:);
    
    % ERROR INDICATORS
    EIs = zeros(1,Ncinit); 
    
    for j = 1:size(Nc0,1)
        disp(j);
        [E,K,L0] = FML_FEM_MATRIX_EXTRACT(Nc0(j,1),[Nc0(j,2),1.24e-03],[Nc0(j,3),1.2e-04],UU);
        ROS = FML_FEM_ROM(E,K,L0,dt,tsim,V);
        ROS_d = (diff(ROS,1,2))./dt;
        ROS_dd = (diff(ROS_d,1,2))./dt;
        ROS_dd = [ROS_dd,zeros(size(K,1),2)];
        RHS = E*ROS_dd + K*ROS;
        Res = L0(2,:) - RHS(2,:);
        NRes = norm(Res);
        Del_T = Tlaplacematrix2d(dt,size(K,1));
        A = (E*Del_T + K);
        EIs(j) = NRes/svds(A,1,'smallest');
    end


    while size(Nc0,1) < N_c
        
        [Imax,Imin,mse,b] = surrogate_FEM(Nc0,EIs);
        Ncn = pointsComputation_FEM(b,Imax,Imin,mse,params);
        Nc0 = [Nc0;Ncn];
        
        for k = Ncinit+1:size(Nc0,1)           
            disp(k)
            [E,K,L0] = FML_FEM_MATRIX_EXTRACT(Nc0(k,1),[Nc0(k,2),1.24e-03],[Nc0(k,3),1.2e-04],UU);
            ROS = FML_FEM_ROM(E,K,L0,dt,tsim,V);
            ROS_d = (diff(ROS,1,2))./dt;
            ROS_dd = (diff(ROS_d,1,2))./dt;
            ROS_dd = [ROS_dd,zeros(size(K,1),2)];
            RHS = E*ROS_dd + K*ROS;
            Res = L0(2,:) - RHS(2,:);
            NRes = norm(Res);
            Del_T = Tlaplacematrix2d(dt,size(K,1));
            A = (E*Del_T + K);
            EIs(k) = NRes/svds(A,1,'smallest'); 
        end
        
        Ncinit = Ncinit + size(Ncn,1);
        
    end

    disp(size(EIs));
    [~,EIsIdx] = max(EIs);
    muNiter = Nc0(EIsIdx,:);
    muNiters(i,:) = muNiter;
    
    % USING SURROGATE MODEL FOR ESTIMATING ERROR
    [Imax,Imin,mse,b] = surrogate_FEM(Nc0,EIs);
    apostError(i) = b(1) + b(2).*muNiter(1) + b(3).*muNiter(2) + b(4).*muNiter(3) + b(5).*muNiter(1)*muNiter(2) + b(6).*muNiter(1)*muNiter(3) + b(7).*muNiter(x2)*muNiter(x3) + b(8).*muNiter(1)*muNiter(2)*muNiter(3) +  b(9).*muNiter(1).^2 + b(10).*muNiter(2).^2 + b(11).*muNiter(3)^2;
    
    vars={'Nc0','Ncinit','A','ET','EIs','HDM','ROS','u','s','v'};
    clear(vars{:});
    
    [E,K,L0] = FML_FEM_MATRIX_EXTRACT(muNiter(1),[muNiter(2),1.24e-03],[muNiter(3),1.2e-04],UU);
    HDM = FML_FEM_HDM(muNiter(1),[muNiter(2),1.24e-03],[muNiter(3),1.2e-04]);
    ROS = FML_FEM_ROM(E,K,L0,dt,tsim,V);
    DiffMat = HDM - ROS;
    
    pointError(i) = norm(HDM(2,:) - ROS(2,:));
    [u_mu,~,~] = svd(DiffMat,'econ');
    
    V(:,i) = u_mu(:,1);
    V = orth(V);
    size(V)
    
end

