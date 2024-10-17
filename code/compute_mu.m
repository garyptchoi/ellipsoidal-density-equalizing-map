function mu = compute_mu(v,f,map)
% Compute the element-wise Beltrami coefficient mu on every triangle
% element. Note that only the norm of mu should be used.
%
% If you use this code in your work, please cite the following paper:
%    Z. Lyu, L. M. Lui, and G. P. T. Choi,
%    "Ellipsoidal Density-Equalizing Map for Genus-0 Closed Surfaces."
%    Preprint, arXiv:2410.12331, 2024.
%
% Copyright (c) 2024, Zhiyuan Lyu, Lok Ming Lui, Gary P. T. Choi
%
% https://github.com/garyptchoi/ellipsoidal-density-equalizing-map

nf = length(f);

% Flatten all triangles of (v,f) onto the plane
v1 = zeros(nf,2);
v2 = [sqrt(sum((v(f(:,1),:) - v(f(:,2),:)).^2,2)),zeros(nf,1)];
A = sum((v(f(:,2),:) - v(f(:,1),:)).*(v(f(:,3),:) - v(f(:,1),:)),2);
B = sqrt(sum((v(f(:,2),:) - v(f(:,1),:)).^2,2));
C = sqrt(sum((v(f(:,3),:) - v(f(:,1),:)).^2,2));
theta_temp = A./(B.*C);
v3 = [sqrt(sum((v(f(:,3),:) - v(f(:,1),:)).^2,2)).*theta_temp,sqrt(sum((v(f(:,3),:) - v(f(:,1),:)).^2,2)).*sqrt(1-theta_temp.^2)];
v_temp = [v1,v2,v3];

% Flatten all triangles of (map,f) onto the plane
map1 = zeros(nf,2);
map2 = [sqrt(sum((map(f(:,1),:) - map(f(:,2),:)).^2,2)),zeros(nf,1)];
A = sum((map(f(:,2),:) - map(f(:,1),:)).*(map(f(:,3),:) - map(f(:,1),:)),2);
B = sqrt(sum((map(f(:,2),:) - map(f(:,1),:)).^2,2));
C = sqrt(sum((map(f(:,3),:) - map(f(:,1),:)).^2,2));
phi_temp = A./(B.*C);
map3 = [sqrt(sum((map(f(:,3),:) - map(f(:,1),:)).^2,2)).*phi_temp,sqrt(sum((map(f(:,3),:) - map(f(:,1),:)).^2,2)).*sqrt(1-phi_temp.^2)];
map_temp = [map1,map2,map3];

% Compute the Beltrami cofficient for every pair of triangles
V = [v_temp(:,1:2);v_temp(:,3:4);v_temp(:,5:6)];   
Map = [map_temp(:,1:2);map_temp(:,3:4);map_temp(:,5:6)]; 
f1 = [(1:nf)',(1:nf)'+nf,(1:nf)'+2*nf];
mu = beltrami_coefficient(V,f1,Map);
