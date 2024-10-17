# Ellipsoidal Density-Equalizing Map

<img src = "https://github.com/garyptchoi/ellipsoidal-density-equalizing-map/blob/main/cover.jpg" height="400" />

* Ellipsoidal density-equalizing map (EDEM): Compute an ellipsoidal density-equalizing map of a genus-0 closed surface onto a prescribed ellipsoid.

* Ellipsoidal density-equalizing quasi-conformal map (EDEQ): Compute an ellipsoidal density-equalizing quasi-conformal map of a genus-0 closed surface, with the ellipsoidal shape automatically determined.

Any comments and suggestions are welcome. 

If you use this code in your work, please cite the following paper:

Z. Lyu,  L. M. Lui, and G. P. T. Choi,
"[Ellipsoidal Density-Equalizing Map for Genus-0 Closed Surfaces.](https://arxiv.org/abs/2410.12331)"
Preprint, arXiv:2410.12331, 2024.

Copyright (c) 2024, Zhiyuan Lyu, Lok Ming Lui, Gary P. T. Choi

https://github.com/garyptchoi/ellipsoidal-density-equalizing-map

===============================================================

EDEM usage:
* `map = EDEM(v,f,population,a,b,c,r0,dt,epsilon,max_iter)`

Input:
* `v`: nv x 3 vertex coordinates of a genus-0 closed surface mesh
* `f`: nf x 3 triangulations of a genus-0 closed surface mesh
* `population`: nf x 1 positive quantity
* `a,b,c`: the radii of the ellipsoid
* `r0`: nv x 3 vertex coordinates of the initial ellipsoidal conformal parameterization 
* `dt`: step size
* `epsilon`: stopping parameter
* `max_iter`: maximum number of iterations

Output:
* `map`: nv x 3 vertex coordinates of the ellipsoidal density-equalizing map

EDEQ usage:
* `[map,a,b,c] = EDEQ(v,f,population,alpha,a0,b0,c0,r0,dt,epsilon,max_iter,shape_iter)`

Input:
* `v`: nv x 3 vertex coordinates of a genus-0 closed surface mesh
* `f`: nf x 3 triangulations of a genus-0 closed surface mesh
* `population`: nf x 1 positive quantity
* `a0,b0,c0`: initial guess of the radii of the ellipsoid
* `r0`: nv x 3 vertex coordinates of the initial ellipsoidal conformal parameterization 
* `dt`: step size
* `epsilon`: stopping parameter
* `max_iter`: maximum number of iterations
* `shape_iter`: number of iterations for each set of fixed radii

Output:
* `map`: nv x 3 vertex coordinates of the ellipsoidal density-equalizing quasi-conformal map
