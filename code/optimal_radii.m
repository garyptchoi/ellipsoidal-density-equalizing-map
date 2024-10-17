function [energy_final,a_final,b_final,c_final,r_final] = optimal_radii(a,b,c,step,v,f,population,r,dt,alpha,interval)
% Optimize the radii of the ellipsoid.
%
% Input:
% v: nv x 3 vertex coordinates of a genus-0 closed surface mesh
% f: nf x 3 triangulations of a genus-0 closed surface mesh
% population: nf x 1 positive quantity
% step: the number of iterations
% a,b,c: the radii of the ellipsoid
% r: nv x 3 vertex coordinates of the last iteration 
% dt: step size (optional, default = 0.1)
% alpha: the parameter of the Beltrami coefficient term
% interval: number of iterations for the interval between two optimized shapes
%
% Output:
% energy_final: the final energy
% a_final,b_final,c_final : the optimized radii
% r_final: nv x 3 vertex coordinates
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
r = (sqrt(1./(r(:,1).^2/a^2 + r(:,2).^2/b^2 + r(:,3).^2/c^2))).*r;

r1 = [r(:,1)/a,r(:,2)/b,r(:,3)/c];

if mod(step,interval) > 0 || step < 2 
    a_final = a;
    b_final = b;
    c_final = c;
    energy_final = energy(a,b,c,v,f,population,r,dt,alpha);
    r_final = r;
    return
end
    
% Stepsize
delta_c = 0.1*(0.9^(step/interval));
delta_b = 0.1*(0.9^(step/interval));

if delta_c < 1e-5
    delta_c = 1e-5;
end

if delta_b < 1e-5
    delta_b = 1e-5;
end

c1 = c - delta_c;
c2 = c + delta_c;

b1 = b - delta_b;
b2 = b + delta_b;

[energy_ori,a,b,c]= energy(a,b,c,v,f,population,r,dt,alpha);

% Varying one parameter
[energy_c1,a,b,c1]= energy(a,b,c1,v,f,population,[r(:,1),r(:,2),r(:,3)/c*c1],dt,alpha);
[energy_c2,a,b,c2]= energy(a,b,c2,v,f,population,[r(:,1),r(:,2),r(:,3)/c*c2],dt,alpha);
[energy_b1,a,b1,c]= energy(a,b1,c,v,f,population,[r(:,1),r(:,2)/b*b1,r(:,3)],dt,alpha);
[energy_b2,a,b2,c]= energy(a,b2,c,v,f,population,[r(:,1),r(:,2)/b*b2,r(:,3)],dt,alpha);

% Varying two parameters
[energy_b1c1,a,b1,c1]= energy(a,b1,c1,v,f,population,[r(:,1),r(:,2)/b*b1,r(:,3)/c*c1],dt,alpha);
[energy_b1c2,a,b1,c2]= energy(a,b1,c2,v,f,population,[r(:,1),r(:,2)/b*b1,r(:,3)/c*c2],dt,alpha);
[energy_b2c1,a,b2,c1]= energy(a,b2,c1,v,f,population,[r(:,1),r(:,2)/b*b2,r(:,3)/c*c1],dt,alpha);
[energy_b2c2,a,b2,c2]= energy(a,b2,c2,v,f,population,[r(:,1),r(:,2)/b*b2,r(:,3)/c*c2],dt,alpha);

% Find the optimal direction
E_all = [energy_ori,a,b,c;
        energy_c1,a,b,c1;
        energy_c2,a,b,c2;
        energy_b1,a,b1,c;
        energy_b2,a,b2,c;
        energy_b1c1,a,b1,c1;
        energy_b1c2,a,b1,c2;
        energy_b2c1,a,b2,c1;
        energy_b2c2,a,b2,c2;];

[~,I] = min(E_all(:,1));

energy_final = E_all(I(1),1);
a_final = E_all(I(1),2);
b_final = E_all(I(1),3);
c_final = E_all(I(1),4);

r_final = [r1(:,1)*a_final,r1(:,2)*b_final,r1(:,3)*c_final]; 

end