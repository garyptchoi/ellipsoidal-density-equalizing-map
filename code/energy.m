function [energy,a,b,c]= energy(a,b,c,v,f,population,r,dt,alpha)
% Compute the energy for given a,b,c.
%
% If you use this code in your work, please cite the following paper:
%    Z. Lyu, L. M. Lui, and G. P. T. Choi,
%    "Ellipsoidal Density-Equalizing Map for Genus-0 Closed Surfaces."
%    Preprint, arXiv:2410.12331, 2024.
%
% Copyright (c) 2024, Zhiyuan Lyu, Lok Ming Lui, Gary P. T. Choi
%
% https://github.com/garyptchoi/ellipsoidal-density-equalizing-map

% normalizing r
r_new = (sqrt(1./(r(:,1).^2/a^2 + r(:,2).^2/b^2 + r(:,3).^2/c^2))).*r;

% re-coupling scheme
rho_f = population./face_area(f,r_new);

rho_v = f2v_area(r_new,f)*rho_f;

% update rho 
L = laplace_beltrami(r_new,f);
A = lumped_mass_matrix(r_new,f);
rho_v_temp = (A+dt*L)\(A*rho_v);

% update density gradient
grad_rho_temp_f = compute_gradient_3D(r_new,f,rho_v_temp);
grad_rho_temp_v = f2v_area(r_new,f)*grad_rho_temp_f;

mu = compute_mu(v,f,r_new);

energy = sum(sum((grad_rho_temp_v/max(rho_v)).^2)) + alpha*sum(abs(mu).^2);

end