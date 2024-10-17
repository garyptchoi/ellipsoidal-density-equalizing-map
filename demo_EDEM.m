% EDEM: Computing ellipsoidal density-equalizing maps of genus-0 closed 
% surfaces using the proposed EDEM method.
%
% Main program:
% map = EDEM(v,f,population,a,b,c,r0,dt,epsilon,max_iter)
% 
% Input:
% v: nv x 3 vertex coordinates of a genus-0 closed surface mesh
% f: nf x 3 triangulations of a genus-0 closed surface mesh
% population: nf x 1 positive quantity
% a,b,c: the radii of the ellipsoid
% r0: nv x 3 vertex coordinates of the initial ellipsoidal conformal parameterization 
%     set r0=[] if you want the algorithm to compute it automatically)
% dt: step size (optional, default = 0.1)
% epsilon: stopping parameter (optional, default = 1e-3)
% max_iter: maximum number of iterations (optional, default = 300)
%
% Output:
% map: nv x 3 vertex coordinates of the ellipsoidal density-equalizing map
%
% If you use this code in your work, please cite the following paper:
%    Z. Lyu, L. M. Lui, and G. P. T. Choi,
%    "Ellipsoidal Density-Equalizing Map for Genus-0 Closed Surfaces."
%    Preprint, arXiv:2410.12331, 2024.
%
% Copyright (c) 2024, Zhiyuan Lyu, Lok Ming Lui, Gary P. T. Choi
%
% https://github.com/garyptchoi/ellipsoidal-density-equalizing-map

addpath('code')
addpath(genpath('data'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 1: Mapping an ellipsoidal surface 

load('ellipsoid_discts.mat'); % Fig. 4a
% load('ellipsoid_cts.mat'); % Fig. 4b
% load('ellipsoid1.mat'); % Fig. 5a
% load('ellipsoid2.mat'); % Fig. 5b
% load('ellipsoid3.mat'); % Fig. 5c
% load('sphere.mat'); % Fig. 5d

max_iter = 300;
dt = 0.1;
epsilon = 1e-3;

% For an input ellipsoidal surface, we can directly run the EDEM algorithm 
% with the initial ellipsoidal parameterization being itself:
r0 = v;
map = EDEM(v,f,population,a,b,c,r0,dt,epsilon,max_iter);

%% Evaluate the density-equalizing property
density1 = population./face_area(f,r0);
density2 = population./face_area(f,map);

var(density1/mean(density1),1)
var(density2/mean(density2),1)

plot_mesh_with_density(v,f,density1);
view([60 10])
title('Input surface');
% plot_mesh_with_density(r0,f,density1);
% view([60 10])
plot_mesh_with_density(map,f,density1);
view([60 10])
title('EDEM result');

% for consistent histogram axes
counts1 = histcounts(density1,0:0.4:10); 
counts2 = histcounts(density2,0:0.4:10); 
max_value = max([counts1, counts2]);

figure;
histogram(density1,0:0.4:10,'FaceColor',[201 0 22]/255,'FaceAlpha',0.5); 
xlim([0 10]);
ylim([0 max_value*1.1]);
xlabel('\rho');
ylabel('Number of faces')
set(gca,'FontSize',15);
set(gca,'linewidth',3);
title('Initial density');

figure;
histogram(density2,0:0.4:10,'FaceColor',[201 0 22]/255,'FaceAlpha',0.5); 
xlim([0 10]);
ylim([0 max_value*1.1]);
xlabel('\rho');
ylabel('Number of faces')
set(gca,'FontSize',15);
set(gca,'linewidth',3);
title('Final density');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 2: Ellipsoidal area-preserving parameterization of genus-0 closed surfaces

load('david.mat'); a = 1; b = 1.2; c = 1.4; % Fig. 6a
% load('buddha.mat'); a = 1; b = 1; c = 1.4; % Fig. 6b
% load('hippocampus.mat'); a = 1; b = 1; c = 2.2; % Fig. 6c
% load('twisted_ball.mat'); a = 1; b = 0.8; c = 1.4; % Fig. 6d


max_iter = 300;
dt = 0.1;
epsilon = 1e-3;

% First compute an initial ellipsoidal conformal parameterization
r0 = ellipsoidal_conformal_map(v,f,a,b,c);

% Set the population as the face area for achieving area-preserving map
population = face_area(f,v);

% Run the EDEM algorithm
map = EDEM(v,f,population,a,b,c,r0,dt,epsilon,max_iter);

%%
plot_mesh(v,f);
title('Input surface');
view([100 10])

plot_mesh(r0,f);
title('Initial parameterization');
view([100 10])

plot_mesh(map,f);
title('EDEM result');
view([100 10])


%% Evaluate the area distortion (for area-preserving maps)

d_area1 = area_distortion(v,f,r0);
d_area2 = area_distortion(v,f,map);

[mean(abs(d_area1)),std(abs(d_area1)),mean(abs(d_area2)),std(abs(d_area2))]

% for consistent histogram axes
counts1 = histcounts(d_area1,-3.5:0.2:3.5); 
counts2 = histcounts(d_area2,-3.5:0.2:3.5); 
max_value = max([counts1, counts2]);

figure;
histogram(d_area1,-3.5:0.2:3.5,'FaceColor',[201 0 22]/255,'FaceAlpha',0.5); 
xlim([-3.5 3.5]);
ylim([0 max_value*1.1]);
xlabel('Logged area ratio');
ylabel('Number of faces')
set(gca,'FontSize',15);
set(gca,'linewidth',3);
title('Initial area distortion');

figure;
histogram(d_area2,-3.5:0.2:3.5,'FaceColor',[201 0 22]/255,'FaceAlpha',0.5); 
xlim([-3.5 3.5]);
ylim([0 max_value*1.1]);
xlabel('Logged area ratio');
ylabel('Number of faces')
set(gca,'FontSize',15);
set(gca,'linewidth',3);
title('Final area distortion');

