clear;
clc;
close all;
format long g;

%% SINEBURST EXCITATION

fex = 120000; T = 1/fex; n = 5; dt = T/20; tsim = 5*n*T;
Vamp = 10e-9; %30;
UU = zeros(1,500);
for i = 1:100
    t = i*dt;
    UU(i) = Vamp*sin(2*pi*fex*t)*sin(pi*fex/n*t)^2;
end

%% PARAMETER SPACE AND DEFINITION

de = 50e6:50e6:50e8;        % Stiffness of the damage
dp = 0.02:0.01:0.11;        % X-position of the damage
ds = 0.002:0.001:0.011;     % Length of the damage

[DE,DP,DS] = ndgrid(de,dp,ds);
paramstable = table(DE(:),DP(:),DS(:));
prior = table2array(paramstable);
clear paramstable;

p = [0, 0.022, 0.25, 0.4, 0.5, 0.8, 1];

%% GENERATE MEASUREMENT DATA FOR THE INTERESTED PARAMETER SET (8e8,0.07,0.01)

Damage_E = 8e8;
Damage_pos = [0.07,1.24e-03];
Damage_size = [0.01,1.2e-04];

tic;
[E,K,L0] = FML_FEM_MATRIX_EXTRACT(Damage_E,Damage_pos,Damage_size,UU);
HDM = FML_FEM_HDM(Damage_E,Damage_pos,Damage_size);
Phi = load('AdaptivePOD_550.mat'); Phi = Phi.V;
ROS = FML_FEM_ROM(E,K,L0,dt,tsim,Phi);
toc;

error = normrnd(0,1e-19,1,497);
measurement = ROS(15648,4:end) + error;
samplings = 7000;
samples = zeros(samplings,3);

%% MCMC - METROPOLIS HASTINGS ALGORITHM

warn = warning('query','all');
id = warn.identifier;
warning('off',id);

theta(1) = de(randi(length(de)));
theta(2) = 0.04;
theta(3) = 0.005;

for ii = 1:samplings

    disp(ii);
    thetaCand(1) = abs(normrnd(theta(1),10e8)); 
    if thetaCand(1)>50e8 
        thetaCand(1) = 10e8 + (50e8 - 10e8).*rand;
    elseif thetaCand(1)<50e6 
        thetaCand(1) = 70e6 + (20e7 - 70e6).*rand;
    end

    thetaCand(2) = round(abs(normrnd(theta(2),0.01)),2);
    if thetaCand(2)>0.11
        thetaCand(2) = round(0.08 + (0.11 - 0.08).*rand,2);
    elseif thetaCand(2)<0.02 
        thetaCand(2) = round(0.03 + (0.06 - 0.03).*rand,2);
    end

    thetaCand(3) = round(abs(normrnd(theta(3),0.001)),3);
    if thetaCand(3)>0.011 
        thetaCand(3) = round(0.008 + (0.011 - 0.008).*rand,3);
    elseif thetaCand(3)<0.002 
        thetaCand(3) = round(0.003 + (0.006 - 0.003).*rand,3);
    end

    if ii == 1
        [E_t,K_t,L0_t] = FML_FEM_MATRIX_EXTRACT(theta(1),[theta(2),1.24e-03],[theta(3),1.2e-04],UU);
        output_theta = FML_FEM_ROM(E_t,K_t,L0_t,dt,tsim,Phi);
        OT = output_theta(15648,4:end);
    end
    
        [E_tc,K_tc,L0_tc] = FML_FEM_MATRIX_EXTRACT(thetaCand(1),[thetaCand(2),1.24e-03],[thetaCand(3),1.2e-04],UU);
        output_thetaCand = FML_FEM_ROM(E_tc,K_tc,L0_tc,dt,tsim,Phi);
        OTC = output_thetaCand(15648,4:end);

    diff = norm((measurement - OT),2);   
    diff_Cand = norm((measurement - OTC),2);

    likelihood = normpdf(diff^0.15,0,0.00027);
    likelihood_Cand = normpdf(diff_Cand^0.15,0,0.00027);

    prior_pdf_de = unifpdf(theta(1),50e6,50e8);
    prior_pdf_dp = unifpdf(theta(2),0.02,0.11);
    prior_pdf_ds = unifpdf(theta(3),0.002,0.011);

    cand_pdf_de = normpdf(thetaCand(1),theta(1),10e8);      
    cand_pdf_dp = normpdf(thetaCand(2),theta(2),0.01);     
    cand_pdf_ds = normpdf(thetaCand(3),theta(3),0.001);     

    posterior = likelihood * prior_pdf_de * prior_pdf_dp * prior_pdf_ds;
    posterior_Cand = likelihood_Cand * cand_pdf_de * cand_pdf_dp * cand_pdf_ds;

    acceptance_ratio = posterior_Cand/posterior;

    if acceptance_ratio > rand(1,1)
        theta(1) = thetaCand(1);
        theta(2) = thetaCand(2); 
        theta(3) = thetaCand(3);
        OT = OTC;
    end
    
    samples(ii,:) = theta;
    toc;
end





