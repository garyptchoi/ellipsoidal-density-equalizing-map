function [map,a,b,c] = EDEQ(v,f,population,alpha,a0,b0,c0,r0,dt,epsilon,max_iter,shape_iter)
% Computing ellipsoidal density-equalizing quasiconformal maps of genus-0 
% closed surfaces using the proposed EDEQ method, with the radii of the
% ellipsoid optimized.
%
% Input:
% v: nv x 3 vertex coordinates of a genus-0 closed surface mesh
% f: nf x 3 triangulations of a genus-0 closed surface mesh
% population: nf x 1 positive quantity
% alpha: the parameter of the Beltrami coefficient term
% a0,b0,c0: initial guess of the radii of the ellipsoid
% r0: nv x 3 vertex coordinates of the initial ellipsoidal conformal parameterization 
%     (set r0=[] if you want the algorithm to compute it automatically)
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

if isempty(r0)
    % Compute an initial ellipsoidal conformal parameterization using the 
    % method in [Choi, J. Comput. Appl. Math. 2024] 
    r0 = ellipsoidal_conformal_map(v,f,a0,b0,c0);
end

if nargin < 9
    dt = 0.1;
end

if nargin < 10
    epsilon = 1e-5;
end

if nargin < 11
    max_iter = 300;
end

if nargin < 12
    shape_iter = 5;
end

r = r0;

% compute density
rho_f = population./face_area(f,r);
rho_v = f2v_area(r,f)*rho_f;

disp('Step    |(E_new-E_old)/E_old|');

bigtri = north_pole_ellipsoid(r0,f);
a = a0;
b = b0;
c = c0;
step = 0;
E_old = energy(a,b,c,r0,f,population,r,dt,alpha);
E_diff = Inf; % for initialization

while E_diff >= epsilon && step < max_iter  
    
    % update rho 
    L = laplace_beltrami(r,f);
    A = lumped_mass_matrix(r,f);
    rho_v_temp = (A+dt*L)\(A*rho_v);

    % update density gradient
    grad_rho_temp_f = compute_gradient_3D(r,f,rho_v_temp);
    grad_rho_temp_v = f2v_area(r,f)*grad_rho_temp_f;
    
    % update displacement
    dr = -[grad_rho_temp_v(:,1)./rho_v_temp,...
           grad_rho_temp_v(:,2)./rho_v_temp,...
           grad_rho_temp_v(:,3)./rho_v_temp];
       
    normal_vector = normal_velocity(r,a,b,c);  
    normal = normal_vector./(sqrt(sum(normal_vector.^2,2)));
    
    % tangential velocity
    dr_proj = dr-[sum(dr.*normal,2),sum(dr.*normal,2),sum(dr.*normal,2)].*normal;
    
    % update and correct overlaps
    r = update_and_correct_overlap_ellipsoid_2(f,r0,r,bigtri,a0,b0,c0,a,b,c,dr_proj,dt);
    
    % update the elliptic radii
    [E_new,a,b,c,r] = optimal_radii(a,b,c,step,r0,f,population,r,dt,alpha,shape_iter);
    
    % evaluate the energy difference
    E_diff = abs(1 - E_new/E_old);
    E_old = E_new;
    step = step + 1;
    disp([num2str(step), '             ',num2str(E_diff)]);

    % re-coupling scheme
    rho_f = population./face_area(f,r);
    rho_v = f2v_area(r,f)*rho_f;
    
end

map = r;

end
