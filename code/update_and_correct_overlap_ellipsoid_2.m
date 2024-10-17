function r_new = update_and_correct_overlap_ellipsoid_2(f,S,r,bigtri,a1,b1,c1,a,b,c,dr,dt)
% Overlap correction scheme for ellipsoidal map with varying elliptic radii.
%
% If you use this code in your work, please cite the following paper:
%    Z. Lyu, L. M. Lui, and G. P. T. Choi,
%    "Ellipsoidal Density-Equalizing Map for Genus-0 Closed Surfaces."
%    Preprint, arXiv:2410.12331, 2024.
%
% Copyright (c) 2024, Zhiyuan Lyu, Lok Ming Lui, Gary P. T. Choi
%
% https://github.com/garyptchoi/ellipsoidal-density-equalizing-map

delta = 0.1; % truncation parameter

r_ori = r;
f_ori = f;
flag = 1;
while flag
    r = r_ori;
    r = r + dt*dr;
    
    % normalization (ensuring x^2/a^2+y^2/b^2+z^2/c^2 = 1)
    r = (sqrt(1./(r(:,1).^2/a^2 + r(:,2).^2/b^2 + r(:,3).^2/c^2))).*r;
    
    %% north pole step
    % rotate the initial and current spheres such that bigtri is the north pole
    S_rotN = rotate_sphere(f,[S(:,1)/a1,S(:,2)/b1,S(:,3)/c1],bigtri); % as we just want to check the bijectivity, we simply rescale the ellipsoid to a sphere here
    [r_rotN,~,RN_inv] = rotate_sphere(f,[r(:,1)/a,r(:,2)/b,r(:,3)/c],bigtri);

    % stereographic projection
    p_S_rotN = stereographic_projection(S_rotN);
    p_r_rotN = stereographic_projection(r_rotN);

    f(bigtri,:) = []; % puncture the triangle

    % ignore the outermost triangles
    [~,I] = sort(S_rotN(:,3),'descend');
    ig_N = I(1:max(round(length(S)/10),3));
    ignore_index_N = find(ismember(f(:,1),ig_N)|ismember(f(:,2),ig_N)|ismember(f(:,3),ig_N));

    % compute Beltrami coefficients
    mu_N = beltrami_coefficient(p_S_rotN,f,p_r_rotN);

    % check if there is any overlap at the inner region
    overlap_N = setdiff(find(abs(mu_N)>=1),ignore_index_N);
    if isempty(overlap_N)
        % no overlap, keep the sphere unchanged
        r_newN = r_rotN;
        north_success = 1;
    else
        % truncation
        mu_N(overlap_N) = (1-delta)*mu_N(overlap_N)./abs(mu_N(overlap_N));

        % LBS
        p_lbsN = linear_beltrami_solver(p_S_rotN,f,mu_N,ig_N,p_r_rotN(ig_N,:));
        mu_N = beltrami_coefficient(p_S_rotN,f,p_lbsN);
        overlap_N = setdiff(find(abs(mu_N)>=1),ignore_index_N);

        if isempty(overlap_N)    
            north_success = 1;
        else
            % still contain overlap, may need a smaller step size
            dt = dt/2;
            north_success = 0;
        end
        r_newN = stereographic_projection(p_lbsN);
    end

    % rotate the sphere to the original position after correction
    r_newN = (RN_inv*r_newN')';

    f = f_ori;
    
    if north_success
        %% south pole step
        % rotate the initial and current spheres such that south_f is the north pole
        south_f = south_pole(f,[r(:,1)/a,r(:,2)/b,r(:,3)/c],bigtri);

        S_rotS = rotate_sphere(f,[S(:,1)/a1,S(:,2)/b1,S(:,3)/c1],south_f);
        [r_rotS,~,RS_inv] = rotate_sphere(f,r_newN,south_f);

        p_S_rotS = stereographic_projection(S_rotS);
        p_r_rotS = stereographic_projection(r_rotS);

        f(south_f,:) = []; % puncture the triangle

        % ignore the outermost triangles
        [~,I1] = sort(S_rotS(:,3),'descend');
        ig_S = I1(1:max(round(length(S)/10),3));
        ignore_index_S = find(ismember(f(:,1),ig_S)|ismember(f(:,2),ig_S)|ismember(f(:,3),ig_S));

        % compute Beltrami coefficients
        mu_S = beltrami_coefficient(p_S_rotS,f,p_r_rotS);

        % check if there is any overlap at the inner region
        overlap_S = setdiff(find(abs(mu_S)>=1),ignore_index_S);
        if isempty(overlap_S)
            % no overlap, keep the sphere unchanged
            r_newS = r_rotS;
            south_success = 1;
        else
            % truncation
            mu_S(overlap_S) = (1-delta)*mu_S(overlap_S)./abs(mu_S(overlap_S));

            % LBS
            p_lbsS = linear_beltrami_solver(p_S_rotS,f,mu_S,ig_S,p_r_rotS(ig_S,:));
            mu_S = beltrami_coefficient(p_S_rotS,f,p_lbsS);
            overlap_S = setdiff(find(abs(mu_S)>=1),ignore_index_S);
            
            if isempty(overlap_S)
                south_success = 1;
            else
                % still contain overlap, may need a smaller step size
                dt = dt/2;
                south_success = 0;
            end
            
            r_newS = stereographic_projection(p_lbsS);
        end
    end
    
    if north_success && south_success
        % done, no need to continue
        flag = 0;
    end
end

% rotate the sphere to the original position after correction
r_new = (RS_inv*r_newS')';
r_new = [r_new(:,1)*a, r_new(:,2)*b, r_new(:,3)*c];



