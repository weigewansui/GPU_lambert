function [V1_1, V1_2, V1_3, V2_1, V2_2, V2_3, extremal_distances_L, extremal_distances_H, exitflag] = ...
        lambert_GPU(r1vec_1, r1vec_2, r1vec_3, r2vec_1, r2vec_2, r2vec_3, tf, m, muC)
%{
LAMBERT_LANCASTERBLANCHARD       High-Thrust Lambert-targeter

lambert_LancasterBlanchard() uses the method developed by 
Lancaster & Blancard, as described in their 1969 paper. Initial values, 
and several details of the procedure, are provided by R.H. Gooding, 
as described in his 1990 paper. 
%}

% Please report bugs and inquiries to: 
%
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com    (personal)
%              oldenhuis@luxspace.lu  (professional)
% Affiliation: LuxSpace s�rl
% Licence    : BSD


% If you find this work useful, please consider a donation:
% https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=6G3S5UYM7HJ3N
      
    % ADJUSTED FOR EML-COMPILATION 29/Sep/2009
    
    % manipulate input
    tol     = 1e-12;                            % optimum for numerical noise v.s. actual precision
    r1      = sqrt(r1vec_1*r1vec_1 + r1vec_2*r1vec_2 + r1vec_3*r1vec_3);              % magnitude of r1vec
    r2      = sqrt(r2vec_1*r2vec_1 + r2vec_2*r2vec_2 + r2vec_3*r2vec_3);              % magnitude of r2vec    
    r1unit_1  = r1vec_1/r1;                         % unit vector of r1vec        
    r1unit_2  = r1vec_2/r1;
    r1unit_3  = r1vec_3/r1;
    
    r2unit_1  = r2vec_1/r2;                         % unit vector of r1vec        
    r2unit_2  = r2vec_2/r2;
    r2unit_3  = r2vec_3/r2;
    
    % unit vector of r2vec
    crsprod_1 = -r1vec_3*r2vec_2 + r1vec_2*r2vec_3;
    crsprod_2 = r1vec_3*r2vec_1 - r1vec_1*r2vec_3;
    crsprod_3 = -r1vec_2*r2vec_1 + r2vec_2*r1vec_1;
    
    mcrsprd = crsprod_1*crsprod_1 + crsprod_2*crsprod_2 + crsprod_3*crsprod_3;   
    
   
    % magnitude of that cross product
    th1unit_1 = -crsprod_3/mcrsprd*r1unit_2 + r1unit_3*crsprod_2/mcrsprd;
    th1unit_2 = -crsprod_1/mcrsprd*r1unit_3 + r1unit_1*crsprod_3/mcrsprd;
    th1unit_3 = -crsprod_3/mcrsprd*r1unit_1 + r1unit_3*crsprod_1/mcrsprd;
    
    th2unit_1 = -crsprod_3/mcrsprd*r2unit_2 + r2unit_3*crsprod_2/mcrsprd;
    th2unit_2 = -crsprod_1/mcrsprd*r2unit_3 + r2unit_1*crsprod_3/mcrsprd;
    th2unit_3 = -crsprod_3/mcrsprd*r2unit_1 + r2unit_3*crsprod_1/mcrsprd;
    
%     th2unit = cross(crsprod/mcrsprd, r2unit);   
    % make 100.4% sure it's in (-1 <= x <= +1)
   
    r1vec_dot_r2vec = r1vec_1*r2vec_1 + r1vec_2*r2vec_2 + r1vec_3*r2vec_3;
    dth = acos( max(-1, min(1, r1vec_dot_r2vec/r1/r2)) ); % turn angle
        
    % if the long way was selected, the turn-angle must be negative
    % to take care of the direction of final velocity
    longway = sign(tf); tf = abs(tf);
    if (longway < 0), dth = dth-2*pi; end
    
    % left-branch
    leftbranch = sign(m); m = abs(m);
    
    % define constants
    c  = sqrt(r1^2 + r2^2 - 2*r1*r2*cos(dth));
    s  = (r1 + r2 + c) / 2;
    T  = sqrt(8*muC/s^3) * tf;
    q  = sqrt(r1*r2)/s * cos(dth/2);
    
    
    % general formulae for the initial values (Gooding)
    % -------------------------------------------------
    
    % some initial values
    
    %%
%     [T, Tp, Tpp, Tppp] = LancasterBlanchard(x, q, m)
%     T0  = LancasterBlanchard(0, q, m);

        % protection against idiotic input
    x = - 2;
    
    % compute parameter E
    E  = x*x - 1;    
    
    % T(x), T'(x), T''(x)
%     if x == 1 % exactly parabolic; solutions known exactly
%         % T(x)
%         T = 4/3*(1-q^3);
        
%     elseif abs(x-1) < 1e-2 % near-parabolic; compute with series
        % evaluate sigma

        
%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          persistent an_1 an_2 ...
%              an_3 an_4 an_5 an_6 an_7 ...
%              an_8 an_9 an_10 an_11 an_12 ...
%              an_13 an_14 an_15 an_16 an_17 ... 
%              an_18 an_19 an_20 an_21 an_22 an_23 an_24 an_25;
            
            an_1 = 4.000000000000000e-001;
            an_2 = 2.142857142857143e-001; 
            an_3 = 4.629629629629630e-002;
            
            an_4 = 6.628787878787879e-003;
            an_5 = 7.211538461538461e-004;
            an_6 = 6.365740740740740e-005;
            an_7 = 4.741479925303455e-006;
            an_8 = 3.059406328320802e-007;
            an_9 = 1.742836409255060e-008;
            an_10 = 8.892477331109578e-010;
            an_11 = 4.110111531986532e-011;
            an_12 = 1.736709384841458e-012;
            an_13 = 6.759767240041426e-014;
            an_14 = 2.439123386614026e-015;
            an_15 = 8.203411614538007e-017;
            an_16 = 2.583771576869575e-018;
            an_17 = 7.652331327976716e-020;
            an_18 = 2.138860629743989e-021;
            an_19 = 5.659959451165552e-023;
            an_20 = 1.422104833817366e-024;
            an_21 = 3.401398483272306e-026;
            an_22 = 7.762544304774155e-028;
            an_23 = 1.693916882090479e-029;
            an_24 = 3.541295006766860e-031;
            an_25 = 7.105336187804402e-033;
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%
    %      expand function 
    %      [sig1, dsigdx1, d2sigdx21, d3sigdx31] = sigmax(-E);
          
            y = -E;

            powers_1 = y^1;
            powers_2 = y^2;
            powers_3 = y^3;
            powers_4 = y^4;
            powers_5 = y^5;
            powers_6 = y^6;
            powers_7 = y^7;
            powers_8 = y^8;
            powers_9 = y^9;
            powers_10 = y^10;
            powers_11 = y^11;
            powers_12 = y^12;
            powers_13 = y^13;
            powers_14 = y^14;
            powers_15 = y^15;
            powers_16 = y^16;
            powers_17 = y^17;
            powers_18 = y^18;
            powers_19 = y^19;
            powers_20 = y^20;
            powers_21 = y^21;
            powers_22 = y^22;
            powers_23 = y^23;
            powers_24 = y^24;
            powers_25 = y^25;

    
             % sigma itself
             sig = 4/3 ...
            + powers_1*an_1 ...
            + powers_2*an_2 ...
            + powers_3*an_3 ...
            + powers_4*an_4 ...
            + powers_5*an_5 ...
            + powers_6*an_6 ...
            + powers_7*an_7 ...
            + powers_8*an_8 ...
            + powers_9*an_9 ...
            + powers_10*an_10 ...
            + powers_11*an_11 ...
            + powers_12*an_12 ...
            + powers_13*an_13 ...
            + powers_14*an_14 ...
            + powers_15*an_15 ...
            + powers_16*an_16 ...
            + powers_17*an_17 ...
            + powers_18*an_18 ...
            + powers_19*an_19 ...
            + powers_20*an_20 ...
            + powers_21*an_21 ...
            + powers_22*an_22 ...
            + powers_23*an_23 ...
            + powers_24*an_24 ...
            + powers_25*an_25; 
    % dsigma / dx (derivative)
%     dsigdx = ( (1:25).*[1, powers(1:24)] ) * an;
      dsigdx = 1 * 1 * an_1 ...
             + 2 * powers_1 * an_2 ...
             + 3 * powers_2 * an_3 ...
             + 4 * powers_3 * an_4 ...
             + 5 * powers_4 * an_5 ...
             + 6 * powers_5 * an_6 ... 
             + 7 * powers_6 * an_7 ...
             + 8 * powers_7 * an_8 ...
             + 9 * powers_8 * an_9 ...
             + 10 * powers_9 * an_10 ...
             + 11 * powers_10 * an_11 ...
             + 12 * powers_11 * an_12 ...
             + 13 * powers_12 * an_13 ...
             + 14 * powers_13 * an_14 ...
             + 15 * powers_14 * an_15 ...
             + 16 * powers_15 * an_16 ...
             + 17 * powers_16 * an_17 ...
             + 18 * powers_17 * an_18 ...
             + 19 * powers_18 * an_19 ...
             + 20 * powers_19 * an_20 ...
             + 21 * powers_20 * an_21 ...
             + 22 * powers_21 * an_22 ...
             + 23 * powers_22 * an_23 ...
             + 24 * powers_23 * an_24 ...
             + 25 * powers_24 * an_25; 
    
    % d2sigma / dx2 (second derivative)
%     d2sigdx2 = ( (1:25).*(0:24).*[1/y, 1, powers(1:23)] ) * an;
    
       d2sigdx2  = 1 * 0 * 1/y * an_1 ...
             + 2 * 1 * 1 * an_2 ...
             + 3 * 2 * powers_1 * an_3 ...
             + 4 * 3 * powers_2 * an_4 ...
             + 5 * 4 * powers_3 * an_5 ...
             + 6 * 5 * powers_4 * an_6 ...
             + 7 * 6 * powers_5 * an_7 ...
             + 8 * 7 * powers_6 * an_8 ...
             + 9 * 8 * powers_7 * an_9 ...
             + 10 * 9 * powers_8 * an_10 ...
             + 11 * 10 * powers_9 * an_11 ...
             + 12 * 11 * powers_10 * an_12 ...
             + 13 * 12 * powers_11 * an_13 ...
             + 14 * 13 * powers_12 * an_14 ...
             + 15 * 14 * powers_13 * an_15 ...
             + 16 * 15 * powers_14 * an_16 ...
             + 17 * 16 * powers_15 * an_17 ...
             + 18 * 17 * powers_16 * an_18 ...
             + 19 * 18 * powers_17 * an_19 ...
             + 20 * 19 * powers_18 * an_20 ...
             + 21 * 20 * powers_19 * an_21 ...
             + 22 * 21 * powers_20 * an_22 ...
             + 23 * 22 * powers_21 * an_23 ...
             + 24 * 23 * powers_22 * an_24 ...
             + 25 * 24 * powers_23 * an_25;
    
    % d3sigma / dx3 (third derivative)
%     d3sigdx3 = ( (1:25).*(0:24).*(-1:23).*[1/y/y, 1/y, 1, powers(1:22)] ) * an;
    
    d3sigdx3 = 1 * 0 * -1 * 1/y/y * an_1 ...
             + 2 * 1 * 0 *  1/y * an_2 ...
             + 3 * 2 * 1 * 1 * an_3 ...
             + 4 * 3 * 2 * powers_1 * an_4 ...
             + 5 * 4 * 3 * powers_2 * an_5 ...
             + 6 * 5 * 4 * powers_3 * an_6 ...
             + 7 * 6 * 5 * powers_4 * an_7 ...
             + 8 * 7 * 6 * powers_5 * an_8 ...
             + 9 * 8 * 7 * powers_6 * an_9 ...
             + 10 * 9 * 8 * powers_7 * an_10 ...
             + 11 * 10 * 9 * powers_8 * an_11 ...
             + 12 * 11 * 10 * powers_9 * an_12 ...
             + 13 * 12 * 11 * powers_10 * an_13 ...
             + 14 * 13 * 12 * powers_11 * an_14 ...
             + 15 * 14 * 13 * powers_12 * an_15 ...
             + 16 * 15 * 14 * powers_13 * an_16 ...
             + 17 * 16 * 15 * powers_14 * an_17 ...
             + 18 * 17 * 16 * powers_15 * an_18 ...
             + 19 * 18 * 17 * powers_16 * an_19 ...
             + 20 * 19 * 18 * powers_17 * an_20 ...
             + 21 * 20 * 19 * powers_18 * an_21 ...
             + 22 * 21 * 20 * powers_19 * an_22 ...
             + 23 * 22 * 21 * powers_20 * an_23 ...
             + 24 * 23 * 22 * powers_21 * an_24 ...
             + 25 * 24 * 23 * powers_22 * an_25;
    
             sig1 = sig;
             dsigdx1 = dsigdx;
             d2sigdx21 = d2sigdx2;
             d3sigdx31 = d3sigdx3;
        %%
%         expand function 
%         [sig2, dsigdx2, d2sigdx22, d3sigdx32] = sigmax(-E*q*q);
        
            y = -E*q*q;

            powers_1 = y^1;
            powers_2 = y^2;
            powers_3 = y^3;
            powers_4 = y^4;
            powers_5 = y^5;
            powers_6 = y^6;
            powers_7 = y^7;
            powers_8 = y^8;
            powers_9 = y^9;
            powers_10 = y^10;
            powers_11 = y^11;
            powers_12 = y^12;
            powers_13 = y^13;
            powers_14 = y^14;
            powers_15 = y^15;
            powers_16 = y^16;
            powers_17 = y^17;
            powers_18 = y^18;
            powers_19 = y^19;
            powers_20 = y^20;
            powers_21 = y^21;
            powers_22 = y^22;
            powers_23 = y^23;
            powers_24 = y^24;
            powers_25 = y^25;

    
             % sigma itself
            sig = 4/3 ...
            + powers_1*an_1 ...
            + powers_2*an_2 ...
            + powers_3*an_3 ...
            + powers_4*an_4 ...
            + powers_5*an_5 ...
            + powers_6*an_6 ...
            + powers_7*an_7 ...
            + powers_8*an_8 ...
            + powers_9*an_9 ...
            + powers_10*an_10 ...
            + powers_11*an_11 ...
            + powers_12*an_12 ...
            + powers_13*an_13 ...
            + powers_14*an_14 ...
            + powers_15*an_15 ...
            + powers_16*an_16 ...
            + powers_17*an_17 ...
            + powers_18*an_18 ...
            + powers_19*an_19 ...
            + powers_20*an_20 ...
            + powers_21*an_21 ...
            + powers_22*an_22 ...
            + powers_23*an_23 ...
            + powers_24*an_24 ...
            + powers_25*an_25; 
    
    % dsigma / dx (derivative)
%     dsigdx = ( (1:25).*[1, powers(1:24)] ) * an;
      dsigdx = 1 * 1 * an_1 ...
             + 2 * powers_1 * an_2 ...
             + 3 * powers_2 * an_3 ...
             + 4 * powers_3 * an_4 ...
             + 5 * powers_4 * an_5 ...
             + 6 * powers_5 * an_6 ... 
             + 7 * powers_6 * an_7 ...
             + 8 * powers_7 * an_8 ...
             + 9 * powers_8 * an_9 ...
             + 10 * powers_9 * an_10 ...
             + 11 * powers_10 * an_11 ...
             + 12 * powers_11 * an_12 ...
             + 13 * powers_12 * an_13 ...
             + 14 * powers_13 * an_14 ...
             + 15 * powers_14 * an_15 ...
             + 16 * powers_15 * an_16 ...
             + 17 * powers_16 * an_17 ...
             + 18 * powers_17 * an_18 ...
             + 19 * powers_18 * an_19 ...
             + 20 * powers_19 * an_20 ...
             + 21 * powers_20 * an_21 ...
             + 22 * powers_21 * an_22 ...
             + 23 * powers_22 * an_23 ...
             + 24 * powers_23 * an_24 ...
             + 25 * powers_24 * an_25; 
    
    % d2sigma / dx2 (second derivative)
%     d2sigdx2 = ( (1:25).*(0:24).*[1/y, 1, powers(1:23)] ) * an;
    
    d2sigdx2  = 1 * 0 * 1/y * an_1 ...
             + 2 * 1 * 1 * an_2 ...
             + 3 * 2 * powers_1 * an_3 ...
             + 4 * 3 * powers_2 * an_4 ...
             + 5 * 4 * powers_3 * an_5 ...
             + 6 * 5 * powers_4 * an_6 ...
             + 7 * 6 * powers_5 * an_7 ...
             + 8 * 7 * powers_6 * an_8 ...
             + 9 * 8 * powers_7 * an_9 ...
             + 10 * 9 * powers_8 * an_10 ...
             + 11 * 10 * powers_9 * an_11 ...
             + 12 * 11 * powers_10 * an_12 ...
             + 13 * 12 * powers_11 * an_13 ...
             + 14 * 13 * powers_12 * an_14 ...
             + 15 * 14 * powers_13 * an_15 ...
             + 16 * 15 * powers_14 * an_16 ...
             + 17 * 16 * powers_15 * an_17 ...
             + 18 * 17 * powers_16 * an_18 ...
             + 19 * 18 * powers_17 * an_19 ...
             + 20 * 19 * powers_18 * an_20 ...
             + 21 * 20 * powers_19 * an_21 ...
             + 22 * 21 * powers_20 * an_22 ...
             + 23 * 22 * powers_21 * an_23 ...
             + 24 * 23 * powers_22 * an_24 ...
             + 25 * 24 * powers_23 * an_25;
             
    % d3sigma / dx3 (third derivative)
%     d3sigdx3 = ( (1:25).*(0:24).*(-1:23).*[1/y/y, 1/y, 1, powers(1:22)] ) * an;
    
        d3sigdx3 = 1 * 0 * -1 * 1/y/y * an_1 ...
             + 2 * 1 * 0 *  1/y * an_2 ...
             + 3 * 2 * 1 * 1 * an_3 ...
             + 4 * 3 * 2 * powers_1 * an_4 ...
             + 5 * 4 * 3 * powers_2 * an_5 ...
             + 6 * 5 * 4 * powers_3 * an_6 ...
             + 7 * 6 * 5 * powers_4 * an_7 ...
             + 8 * 7 * 6 * powers_5 * an_8 ...
             + 9 * 8 * 7 * powers_6 * an_9 ...
             + 10 * 9 * 8 * powers_7 * an_10 ...
             + 11 * 10 * 9 * powers_8 * an_11 ...
             + 12 * 11 * 10 * powers_9 * an_12 ...
             + 13 * 12 * 11 * powers_10 * an_13 ...
             + 14 * 13 * 12 * powers_11 * an_14 ...
             + 15 * 14 * 13 * powers_12 * an_15 ...
             + 16 * 15 * 14 * powers_13 * an_16 ...
             + 17 * 16 * 15 * powers_14 * an_17 ...
             + 18 * 17 * 16 * powers_15 * an_18 ...
             + 19 * 18 * 17 * powers_16 * an_19 ...
             + 20 * 19 * 18 * powers_17 * an_20 ...
             + 21 * 20 * 19 * powers_18 * an_21 ...
             + 22 * 21 * 20 * powers_19 * an_22 ...
             + 23 * 22 * 21 * powers_20 * an_23 ...
             + 24 * 23 * 22 * powers_21 * an_24 ...
             + 25 * 24 * 23 * powers_22 * an_25;
             
             sig2 = sig;
             dsigdx2 = dsigdx;
             d2sigdx22 = d2sigdx2;
             d3sigdx32 = d3sigdx3;
        %%
        
        % T(x)
        T0 = sig1 - q^3*sig2;
        
%     else % all other cases
        % compute all substitution functions
        y  = sqrt(abs(E));
        z  = sqrt(1 + q^2*E);
        f  = y*(z - q*x);
        g  = x*z - q*E;

        % BUGFIX: (Simon Tardivel) this line is incorrect for E==0 and f+g==0
        % d  = (E < 0)*(atan2(f, g) + pi*m) + (E > 0)*log( max(0, f + g) );
        % it should be written out like so: 
        if (E<0)
            d = atan2(f, g) + pi*m;
        elseif (E==0)
            d = 0;
        else 
            d = log(max(0, f+g));
        end

        % T(x)
        T0 = 2*(x - q*z - d/y)/E;
%     end
    %%
    
    

    Td  = T0 - T;
    phr = mod(2*atan2(1 - q^2, 2*q), 2*pi);
    
    % initial output is pessimistic
    V1_1 = NaN;
    V1_2 = NaN;
    V1_3 = NaN;
    
    V2_1 = V1_1;
    V2_2 = V1_2;
    V2_3 = V1_3;
    
    extremal_distances_L = NaN;
    extremal_distances_H = NaN;
    
    % single-revolution case
    if (m == 0)
        x01 = T0*Td/4/T;
        if (Td > 0)
            x0 = x01;
        else
            x01 = Td/(4 - Td);
            x02 = -sqrt( -Td/(T+T0/2) );
            W   = x01 + 1.7*sqrt(2 - phr/pi);
            if (W >= 0)
                x03 = x01;
            else
                x03 = x01 + (-W).^(1/16).*(x02 - x01);
            end
            lambda = 1 + x03*(1 + x01)/2 - 0.03*x03^2*sqrt(1 + x01);
            x0 = lambda*x03;
        end
        
        % this estimate might not give a solution
        if (x0 < -1), exitflag = -1; return; end
        
    % multi-revolution case
    else
        
        % determine minimum Tp(x)
        xMpi = 4/(3*pi*(2*m + 1));        
        if (phr < pi)
            xM0 = xMpi*(phr/pi)^(1/8);
        elseif (phr > pi)
            xM0 = xMpi*(2 - (2 - phr/pi)^(1/8));
        % EMLMEX requires this one
        else
            xM0 = 0;
        end
        
        % use Halley's method
        xM = xM0;  Tp = inf;  iterations = 0;
        while abs(Tp) > tol            
            % iterations
            iterations = iterations + 1;            
            % compute first three derivatives
            %%%%%%%%%%%%%%%%%%%%%%%
%             [dummy, Tp, Tpp, Tppp] = LancasterBlanchard(xM, q, m);%#ok
            
            %%
            %   expand LancasterBlanchard
            
            x = xM;
            
            if (x < -1) % impossible; negative eccentricity
                x = abs(x) - 2;
            elseif (x == -1) % impossible; offset x slightly
                x = x + eps;
            end

            % compute parameter E
            E  = x*x - 1;    

            % T(x), T'(x), T''(x)
            if x == 1 % exactly parabolic; solutions known exactly
                % T(x)
                T = 4/3*(1-q^3);
                % T'(x)
                Tp = 4/5*(q^5 - 1);
                % T''(x)
                Tpp = Tp + 120/70*(1 - q^7);
                % T'''(x)
                Tppp = 3*(Tpp - Tp) + 2400/1080*(q^9 - 1);

            elseif abs(x-1) < 1e-2 % near-parabolic; compute with series
                % evaluate sigma
%                 [sig1, dsigdx1, d2sigdx21, d3sigdx31] = sigmax(-E);
                %%
                y = -E;

                sigm;
    
                sig1 = sig;
                dsigdx1 = dsigdx;
                d2sigdx21 = d2sigdx2;
                d3sigdx31 = d3sigdx3;
             %%
             
%                 [sig2, dsigdx2, d2sigdx22, d3sigdx32] = sigmax(-E*q*q);
                y = -E*q*q;
                sigm;
                sig2 = sig;
                dsigdx2 = dsigdx;
                d2sigdx22 = d2sigdx2;
                d3sigdx32 = d3sigdx3;
                
                % T(x)
                T = sig1 - q^3*sig2;
                % T'(x)
                Tp = 2*x*(q^5*dsigdx2 - dsigdx1);
                % T''(x)        
                Tpp = Tp/x + 4*x^2*(d2sigdx21 - q^7*d2sigdx22);
                % T'''(x)
                Tppp = 3*(Tpp-Tp/x)/x + 8*x*x*(q^9*d3sigdx32 - d3sigdx31);

            else % all other cases
                % compute all substitution functions
                y  = sqrt(abs(E));
                z  = sqrt(1 + q^2*E);
                f  = y*(z - q*x);
                g  = x*z - q*E;

                % BUGFIX: (Simon Tardivel) this line is incorrect for E==0 and f+g==0
                % d  = (E < 0)*(atan2(f, g) + pi*m) + (E > 0)*log( max(0, f + g) );
                % it should be written out like so: 
                if (E<0)
                    d = atan2(f, g) + pi*m;
                elseif (E==0)
                    d = 0;
                else 
                    d = log(max(0, f+g));
                end

                % T(x)
                T = 2*(x - q*z - d/y)/E;
                %  T'(x)
                Tp = (4 - 4*q^3*x/z - 3*x*T)/E;
                % T''(x)
                Tpp = (-4*q^3/z * (1 - q^2*x^2/z^2) - 3*T - 3*x*Tp)/E;
                % T'''(x) 
                Tppp = (4*q^3/z^2*((1 - q^2*x^2/z^2) + 2*q^2*x/z^2*(z - x)) - 8*Tp - 7*x*Tpp)/E;     

            end
            
            dummy = T;
              
            %%
            %%%%%%%%%%%%%%%%%%%%%%%
            % new value of xM
            xMp = xM;
            xM  = xM - 2*Tp*Tpp / (2*Tpp^2 - Tp*Tppp);            
            % escape clause
            if mod(iterations, 7), xM = (xMp+xM)/2; end
            % the method might fail. Exit in that case
            if (iterations > 25), exitflag = -2; return; end
        end
        
        % xM should be elliptic (-1 < x < 1)
        % (this should be impossible to go wrong)
        if (xM < -1) || (xM > 1), exitflag = -1; return; end
        
        % corresponding time
%         TM = LancasterBlanchard(xM, q, m);
        x = xM;
        LancasterBlanchard;
        TM = T;
        
        % T should lie above the minimum T
        if (TM > T), exitflag = -1; return; end
        
        % find two initial values for second solution (again with lambda-type patch)
        % --------------------------------------------------------------------------
        
        % some initial values
        TmTM  = T - TM;   T0mTM = T0 - TM;
%         [dummy, Tp, Tpp] = LancasterBlanchard(xM, q, m);%#ok
        x = xM;
        LancasterBlanchard;
        dummy = T;
                
        % first estimate (only if m > 0)
        if leftbranch > 0
            x   = sqrt( TmTM / (Tpp/2 + TmTM/(1-xM)^2) );
            W   = xM + x;
            W   = 4*W/(4 + TmTM) + (1 - W)^2;
            x0  = x*(1 - (1 + m + (dth - 1/2)) / ...
                (1 + 0.15*m)*x*(W/2 + 0.03*x*sqrt(W))) + xM;
            
            % first estimate might not be able to yield possible solution
            if (x0 > 1), exitflag = -1; return; end
            
        % second estimate (only if m > 0)
        else
            if (Td > 0)
                x0 = xM - sqrt(TM/(Tpp/2 - TmTM*(Tpp/2/T0mTM - 1/xM^2)));
            else
                x00 = Td / (4 - Td);
                W = x00 + 1.7*sqrt(2*(1 - phr));
                if (W >= 0)
                    x03 = x00;
                else
                    x03 = x00 - sqrt((-W)^(1/8))*(x00 + sqrt(-Td/(1.5*T0 - Td)));
                end
                W      = 4/(4 - Td);
                lambda = (1 + (1 + m + 0.24*(dth - 1/2)) / ...
                    (1 + 0.15*m)*x03*(W/2 - 0.03*x03*sqrt(W)));
                x0     = x03*lambda;
            end
            
            % estimate might not give solutions
            if (x0 < -1), exitflag = -1; return; end
            
        end
    end
    
    % find root of Lancaster & Blancard's function
    % --------------------------------------------
    
    % (Halley's method)    
    x = x0;
    Tx = inf; 
    iterations = 0;
    while abs(Tx) > tol        
        % iterations
        iterations = iterations + 1;        
        % compute function value, and first two derivatives
%         [Tx, Tp, Tpp] = LancasterBlanchard(x, q, m); 
        LancasterBlanchard;
        Tx = T;
        
        % find the root of the *difference* between the
        % function value [T_x] and the required time [T]
        Tx = Tx - T;
        % new value of x
        xp = x;
        x  = x - 2*Tx*Tp / (2*Tp^2 - Tx*Tpp);                 
        % escape clause
        if mod(iterations, 7), x = (xp+x)/2; end        
        % Halley's method might fail
        if iterations > 25, exitflag = -2; return; end
    end   
    
    % calculate terminal velocities
    % -----------------------------
    
    % constants required for this calculation
    gamma = sqrt(muC*s/2);
    if (c == 0)
        sigma = 1;
        rho   = 0;
        z     = abs(x);
    else
        sigma = 2*sqrt(r1*r2/(c^2)) * sin(dth/2);
        rho   = (r1 - r2)/c;
        z     = sqrt(1 + q^2*(x^2 - 1));
    end
    
    % radial component
    Vr1    = +gamma*((q*z - x) - rho*(q*z + x)) / r1;
    
    Vr1vec_1 = Vr1*r1unit_1;
    Vr1vec_2 = Vr1*r1unit_2;
    Vr1vec_3 = Vr1*r1unit_3;
    
    Vr2    = -gamma*((q*z - x) + rho*(q*z + x)) / r2;
    
    Vr2vec_1 = Vr2*r2unit_1;
    Vr2vec_2 = Vr2*r2unit_2;
    Vr2vec_3 = Vr2*r2unit_3;
    
    % tangential component
    Vtan1    = sigma * gamma * (z + q*x) / r1;
    Vtan1vec_1 = Vtan1 * th1unit_1;
    Vtan1vec_2 = Vtan1 * th1unit_2;
    Vtan1vec_3 = Vtan1 * th1unit_3;
    
    Vtan2    = sigma * gamma * (z + q*x) / r2;
    Vtan2vec_1 = Vtan2 * th2unit_1;
    Vtan2vec_2 = Vtan2 * th2unit_2;
    Vtan2vec_3 = Vtan2 * th2unit_3;
    
    % Cartesian velocity
    V1_1 = Vtan1vec_1 + Vr1vec_1;
    V1_2 = Vtan1vec_2 + Vr1vec_2;
    V1_3 = Vtan1vec_3 + Vr1vec_3;
    
    V2_1 = Vtan2vec_1 + Vr2vec_1;
    V2_2 = Vtan2vec_2 + Vr2vec_2;
    V2_3 = Vtan2vec_3 + Vr2vec_3;
    
    % exitflag
    exitflag = 1; % (success)
    
    % also determine minimum/maximum distance
    a = s/2/(1 - x^2); % semi-major axis
    
    minmax_distances;
%     
%     minmax_distances(r1vec, r1, r2vec, r2, dth, a, V1, V2, m, muC)
%     r1vec_1
%     = minmax_distances(r1vec, r1, r1vec, r2, dth, a, V1, V2, m, muC);
    
end

% Lancaster & Blanchard's function, and three derivatives thereof


% compute minimum and maximum distances to the central body

