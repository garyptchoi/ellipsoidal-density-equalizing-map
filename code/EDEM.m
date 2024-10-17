function map = EDEM(v,f,population,a,b,c,r0,dt,epsilon,max_iter)
% Computing ellipsoidal density-equalizing maps of genus-0 closed surfaces 
% using the proposed EDEM method.
%
% Input:
% v: nv x 3 vertex coordinates of a genus-0 closed surface mesh
% f: nf x 3 triangulations of a genus-0 closed surface mesh
% population: nf x 1 positive quantity
% a,b,c: the radii of the ellipsoid
% r0: nv x 3 vertex coordinates of the initial ellipsoidal conformal parameterization 
%     (set r0=[] if you want the algorithm to compute it automatically)
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

if isempty(r0)
    % Compute an initial ellipsoidal conformal parameterization using the 
    % method in [Choi, J. Comput. Appl. Math. 2024] 
    r0 = ellipsoidal_conformal_map(v,f,a,b,c);
end

if nargin < 8
    dt = 0.1;
end

if nargin < 9
    epsilon = 1e-3;
end

if nargin < 10
    max_iter = 300;
end

r = r0;
% compute density
rho_f = population./face_area(f,r);
rho_v = f2v_area(r,f)*rho_f;

step = 0;
rho_v_error = std(rho_v)/mean(rho_v);
disp('Step     std(rho)/mean(rho)');
disp([num2str(step), '        ',num2str(rho_v_error)]);

    
bigtri = north_pole_ellipsoid(r0,f);
    
while rho_v_error >= epsilon && step < max_iter
    
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
    r = update_and_correct_overlap_ellipsoid(f,r0,r,bigtri,a,b,c,dr_proj,dt);
    
    step = step + 1;
    rho_v_error = std(rho_v_temp)/mean(rho_v_temp);
    disp([num2str(step), '        ',num2str(rho_v_error)]);

    % re-coupling scheme
    rho_f = population./face_area(f,r);
    rho_v = f2v_area(r,f)*rho_f;
    
end

map = r;

end
