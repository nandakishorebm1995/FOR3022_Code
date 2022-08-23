function [U_COM] = FML_FEM_HDM(Damage_E,Damage_pos,Damage_size)

% This is the forward FEM function that extracts the system matrices
% from COMSOL and solves the second order system of equations using
% Newmark's beta time integration method. 

% COMSOL solves for the actual displacement field itself rather than 
% solving for DOFs, a usual approach in classical FEM.
% The load matrix was generated manually by applying sineburst function at
% the nodes where the excitation is applied. 

% Currently, besides the damage position parameters, only the elasticity
% modulus of the defect is considered for the parameter reconstruction.

%% LOADING MODEL
model = mphload('Stahl_laminat_032021_mesh');
% model.param.set('n', '5', 'Anzahl der Anregungen');
% mphgetexpressions(model.param)

%% DAMAGE PARAMETER DEFINITION

% DAMAGE MATERIAL DEFINITION
model.param.set('DE', strcat(num2str(Damage_E),'[Pa]'), 'Damage youngs modulus');
model.component('comp1').material('mat12').propertyGroup('Enu').set('youngsmodulus', 'DE');

% DAMAGE GEOMETRY DEFINITION
model.component('comp1').geom('geom1').feature('r2').set('pos', {num2str(Damage_pos(1)), num2str(Damage_pos(2))});
model.component('comp1').geom('geom1').feature('r2').set('size', {num2str(Damage_size(1)), num2str(Damage_size(2))});
model.component('comp1').geom('geom1').feature('ls1').set('coord1', [Damage_pos(1), 0]);
model.component('comp1').geom('geom1').feature('ls1').set('coord2', [Damage_pos(1), 0.00198]);
model.component('comp1').geom('geom1').feature('ls2').set('coord1', [Damage_pos(1)+Damage_size(1), 0]);
model.component('comp1').geom('geom1').feature('ls2').set('coord2', [Damage_pos(1)+Damage_size(1), 0.00198]);
model.component('comp1').geom('geom1').run;

%% MESHING

model.component('comp1').mesh('mesh1').feature('size').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('size').set('hmax', 0.0005);
model.component('comp1').mesh('mesh1').feature('size').set('hmin', 0.0005);
model.component('comp1').mesh('mesh1').feature('map1').active(false);
model.component('comp1').mesh('mesh1').feature('map1').feature('dis1').set('numelem', 10);
model.component('comp1').mesh('mesh1').feature('map1').feature('dis2').set('numelem', 2);
model.component('comp1').mesh('mesh1').feature('map1').feature('dis3').set('numelem', 100);
model.component('comp1').mesh('mesh1').run;
model.component('comp1').mesh('mesh1').run;
model.component('comp1').mesh('mesh1').feature('map2').selection.set([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48]);
model.component('comp1').mesh('mesh1').run;
model.sol('sol1').runAll;

U_COM = mphgetu(model,'solnum',[1:1:500]);

%%
% %% EXTRACTION OF SYSTEM MATRICES
% 
% str = mphmatrix(model,'sol1','out',{'E','K','L'},'initmethod','init','initsol','zero');
% E = str.E; 
% K = str.K;
% sp = size(K,1);
% 
% %% EVALUATION OF LOAD MATRIX
% 
% L0 = zeros(sp,500);
% L0(2,:) = UU;
% L0(966,:) = UU;

%% SWITCH THIS ON FOR COMPARISON OF COMSOL AND SUBROUTINE SOLUTION
%% NUMERICAL INTEGRATION OF THE SECOND ORDER ODE (NEWMARK'S BETA METHOD)

% tic;
% [dircoeff] = NewmarkIntegration3(E,K,L0,dt,sp,tsim);
% toc;

% PLOTTING
% Comparison with the solution obtained from COMSOL

% figure;plot(linspace(0,tsim,500),dircoeff(3504,:),'linewidth',2);
% fid = fopen('PlotCompare.txt','rt');
% CD = textscan(fid,'%f%f','Delimiter','\t');
% ComData = [CD{1,1},CD{1,2}];
% 
% hold on; plot(ComData(:,1),ComData(:,2),'linewidth',2);
% grid on; legend('Matlab subroutine solution', 'COMSOL solution');
% title('Wave propagation recorded at (0.0575,0.0014225)');
% xlabel('Time in s');
% ylabel('Vertical displacement in m');

end

