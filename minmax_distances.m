% function extremal_distances = ...
%         minmax_distances(r1vec, r1, r2vec, r2, dth, a, V1, V2, m, muC)
        
    % default - minimum/maximum of r1,r2
    minimum_distance = min(r1_tmp,r2_tmp);
    maximum_distance = max(r1_tmp,r2_tmp);
    
    % was the longway used or not?
    longway = abs(dth_tmp) > pi;
    
    % eccentricity vector (use triple product identity)
    norm_V1 = (V1_1_tmp*V1_1_tmp + V1_2_tmp*V1_2_tmp + V1_3_tmp*V1_3_tmp);
    V1_dot_r1vec = V1_1*r1vec_1 + V1_2*r1vec_2 + V1_3*r1vec_3;
    
    evec_1 = (norm_V1 * r1vec_1_tmp - V1_dot_r1vec*V1_1_tmp)/muC_tmp - r1vec_1_tmp/r1_tmp;

    evec_2 = (norm_V1 * r1vec_2_tmp - V1_dot_r1vec*V1_2_tmp)/muC_tmp - r1vec_2_tmp/r1_tmp;
    evec_3 = (norm_V1 * r1vec_3_tmp - V1_dot_r1vec*V1_3_tmp)/muC_tmp - r1vec_2_tmp/r1_tmp;
    
    % eccentricity
    e = sqrt(evec_1*evec_1 + evec_2*evec_2 + evec_3*evec_3);     
    % apses
    pericenter = a_tmp*(1-e);
    apocenter  = inf;                    % parabolic/hyperbolic case
    if (e < 1), apocenter = a_tmp*(1+e); end % elliptic case
    
    % since we have the eccentricity vector, we know exactly where the
    % pericenter lies. Use this fact, and the given value of [dth], to 
    % cross-check if the trajectory goes past it
    if (m > 0) % obvious case (always elliptical and both apses are traversed)
        minimum_distance = pericenter;
        maximum_distance = apocenter;
    else % less obvious case
        % compute theta1&2 ( use (AxB)-(CxD) = (C·B)(D·A) - (C·A)(B·D) ))
        
        evec_dot_V1 = evec_1*V1_1_tmp + evec_2*V1_2_tmp + evec_3*V1_3_tmp;
        r1vec_dot_evec = r1vec_1_tmp*evec_1 + r1vec_2_tmp*evec_2 +r1vec_3_tmp*evec_3;
        
        pm1 = sign( r1_tmp*r1_tmp*evec_dot_V1 - r1vec_dot_evec* V1_dot_r1vec );  
        
        evec_dot_V2 = evec_1*V2_1_tmp + evec_2*V2_2_tmp + evec_3*V2_3_tmp;
        r2vec_dot_evec = r2vec_1_tmp*evec_1 + r2vec_2_tmp*evec_2 + r2vec_3_tmp*evec_3;
        r2vec_dot_V2 = r2vec_1_tmp*V2_1_tmp + r2vec_2_tmp*V2_2_tmp + r2vec_3_tmp * V2_3_tmp; 
        
        pm2 = sign( r2_tmp*r2_tmp*(evec_dot_V2) - (r2vec_dot_evec )*(r2vec_dot_V2) );  
        % make 100.4% sure it's in (-1 <= theta12 <= +1)

        theta1 = pm1*acos( max(-1, min(1, r1vec_dot_evec/r1/e)) );
        
        theta2 = pm2*acos( max(-1, min(1, r2vec_dot_evec/r2/e)) );
        % points 1&2 are on opposite sides of the symmetry axis -- minimum 
        % and maximum distance depends both on the value of [dth], and both 
        % [theta1] and [theta2]
        if (theta1*theta2 < 0)
            % if |th1| + |th2| = turnangle, we know that the pericenter was 
            % passed
            if (abs(theta1)+abs(theta2) == dth_tmp)
                minimum_distance = pericenter;
            % this condition can only be false for elliptic cases, and 
            % when it is indeed false, we know that the orbit passed 
            % apocenter
            else
                maximum_distance = apocenter;
            end
        % points 1&2 are on the same side of the symmetry axis. Only if the 
        % long-way was used are the min. and max. distances different from 
        % the min. and max. values of the radii (namely, equal to the apses)
        elseif longway
            minimum_distance = pericenter;
            if (e < 1), maximum_distance = apocenter; end
        end
    end
    
    % output argument
    extremal_distances_L_tmp = minimum_distance;
    extremal_distances_H_tmp = maximum_distance;
    