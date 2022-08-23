function [E,K,L0] = FML_FEM_MATRIX_EXTRACT(Damage_E,Damage_pos,Damage_size,UU)

% This function extracts the system matrices from COMSOL

%% LOADING MODEL

model = mphload('Stahl_laminat_032021_mesh');
% figure; mphgeom(model);

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
% model.sol('sol1').runAll;

%% MESHING

model.component('comp1').mesh('mesh1').feature('size').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('size').set('hmax', 0.0005);
model.component('comp1').mesh('mesh1').feature('size').set('hmin', 0.0005);
model.component('comp1').mesh('mesh1').feature('map1').active(false);
model.component('comp1').mesh('mesh1').feature('map1').feature('dis1').set('numelem', 10);
model.component('comp1').mesh('mesh1').feature('map1').feature('dis2').set('numelem', 2);
model.component('comp1').mesh('mesh1').feature('map1').feature('dis3').set('numelem', 100);
model.component('comp1').mesh('mesh1').run;
model.component('comp1').mesh('mesh1').feature('map2').selection.set([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48]);
model.component('comp1').mesh('mesh1').run;
% figure; mphmesh(model);

%% EXTRACTION OF SYSTEM MATRICES

str = mphmatrix(model,'sol1','out',{'E','K','L'},'initmethod','init','initsol','zero');
E = str.E; 
K = str.K;
sp = size(K,1);

%% EVALUATION OF LOAD MATRIX

L0 = zeros(sp,500);
L0(2,:) = UU;
L0(966,:) = UU;

end


