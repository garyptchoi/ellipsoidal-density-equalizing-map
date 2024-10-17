% EDEQ: Computing ellipsoidal density-equalizing quasiconformal maps of 
% genus-0 closed surfaces using the proposed EDEQ method, with the radii of 
% the ellipsoid optimized.

% Main program:
% [map,a,b,c] = EDEQ(v,f,population,alpha,a0,b0,c0,r0,dt,epsilon,max_iter,shape_iter)
% 
% Input:
% v: nv x 3 vertex coordinates of a genus-0 closed surface mesh
% f: nf x 3 triangulations of a genus-0 closed surface mesh
% population: nf x 1 positive quantity
% alpha: the parameter of the Beltrami coefficient term
% a0,b0,c0: initial guess of the radii of the ellipsoid
% r0: nv x 3 vertex coordinates of the initial ellipsoidal conformal parameterization 
%     set r0=[] if you want the algorithm to compute it automatically)
% dt: step size (optional, default = 0.1)
% epsilon: stopping parameter (optional, default = 1e-5)
% max_iter: maximum number of iterations (optional, default = 300)
% shape_iter: number of iterations for each set of fixed radii (optional, default = 5)
%
% Output:
% map: nv x 3 vertex coordinates of the ellipsoidal density-equalizing quasiconformal map
% a,b,c: the optimal radii of the ellipsoid
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

load('ellipsoid_discts_edeq.mat'); % Fig. 7a
% load('ellipsoid_cts_edeq.mat'); % Fig. 7b
% load('ellipsoid1_edeq.mat'); % Fig. 8a
% load('ellipsoid2_edeq.mat'); % Fig. 8b
% load('ellipsoid3_edeq.mat'); % Fig. 8c
% load("sphere_edeq.mat"); % Fig. 8d


alpha = 5;
dt = 0.1;
epsilon = 1e-5;
max_iter = 300;
shape_iter = 5;

% For an input ellipsoidal surface, we can directly run the EDEQ algorithm 
% with the initial ellipsoidal parameterization being itself:
r0 = v;
[map,a,b,c] = EDEQ(v,f,population,alpha,a0,b0,c0,r0,dt,epsilon,max_iter,shape_iter);

%% Evaluate the performance

% Evaluate the density-equalizing property
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
title('EDEQ result');


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

% Evaluate the quasiconformal distortion
mu = compute_mu(r0,f,map);
mean(abs(mu))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 2: Ellipsoidal area-preserving parameterization of genus-0 closed surfaces

% Starting with an initial ellipsoid (Fig. 9):
% load('duck.mat'); a0 = 1; b0 = 1; c0 = 1.2; % Fig. 9a
% load('hippocampus.mat'); a0 = 1; b0 = 1; c0 = 1.5; % Fig. 9b
% load('vaselion.mat'); a0 = 1; b0 = 1; c0 = 1.4; % Fig. 9c

% Starting with an initial sphere (Fig. 10):
load('duck.mat'); a0 = 1; b0 = 1; c0 = 1; % Fig. 10a
% load('hippocampus.mat'); a0 = 1; b0 = 1; c0 = 1; % Fig. 10b
% load('vaselion.mat'); a0 = 1; b0 = 1; c0 = 1; % Fig. 10c

alpha = 5;
dt = 0.1;
epsilon = 1e-5;
max_iter = 300;
shape_iter = 5;

% First compute an initial ellipsoidal conformal parameterization
r0 = ellipsoidal_conformal_map(v,f,a0,b0,c0);

% Set the population as the face area for achieving area-preserving map
population = face_area(f,v);

% Run the EDEQ algorithm
[map,a,b,c] = EDEQ(v,f,population,alpha,a0,b0,c0,r0,dt,epsilon,max_iter,shape_iter);

%%

plot_mesh(v,f);
title('Input surface');
view([100 10])

plot_mesh(r0,f);
title('Initial parameterization');
view([100 10])

plot_mesh(map,f);
title('EDEQ result');
view([100 10])


%% Evaluate the performance

% Evaluate the area distortion (for area-preserving maps)

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

% Evaluate the quasiconformal distortion
mu = compute_mu(v,f,map);
mean(abs(mu))
