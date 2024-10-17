function map = north_pole_ellipsoid(v,f)
% Compute the north pole of the ellipsoid.
%
% If you use this code in your work, please cite the following paper:
%    Z. Lyu, L. M. Lui, and G. P. T. Choi,
%    "Ellipsoidal Density-Equalizing Map for Genus-0 Closed Surfaces."
%    Preprint, arXiv:2410.12331, 2024.
%
% Copyright (c) 2024, Zhiyuan Lyu, Lok Ming Lui, Gary P. T. Choi
%
% https://github.com/garyptchoi/ellipsoidal-density-equalizing-map

f_center = (v(f(:,1),:)+v(f(:,2),:)+v(f(:,3),:))/3;

[~,index] = max(f_center(:,3));
map = index;
end

