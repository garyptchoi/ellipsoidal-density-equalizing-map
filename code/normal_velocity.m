function map = normal_velocity(v,a,b,c)
% Compute the normal vector of the ellipsoid.
%
% If you use this code in your work, please cite the following paper:
%    Z. Lyu, L. M. Lui, and G. P. T. Choi,
%    "Ellipsoidal Density-Equalizing Map for Genus-0 Closed Surfaces."
%    Preprint, arXiv:2410.12331, 2024.
%
% Copyright (c) 2024, Zhiyuan Lyu, Lok Ming Lui, Gary P. T. Choi
%
% https://github.com/garyptchoi/ellipsoidal-density-equalizing-map

map = [(2*v(:,1))/a^2,(2*v(:,2))/b^2,(2*v(:,3))/c^2];
end

