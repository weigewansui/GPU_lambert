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
    
    mcrsprd = sqrt(crsprod_1*crsprod_1 + crsprod_2*crsprod_2 + crsprod_3*crsprod_3);   
    
    % magnitude of that cross product
    th1unit_1 = -crsprod_3/mcrsprd*r1unit_2 + r1unit_3*crsprod_2/mcrsprd;
    th1unit_2 = -crsprod_1/mcrsprd*r1unit_3 + r1unit_1*crsprod_3/mcrsprd;
    th1unit_3 = -crsprod_2/mcrsprd*r1unit_1 + r1unit_2*crsprod_1/mcrsprd;
    
    th2unit_1 = -crsprod_3/mcrsprd*r2unit_2 + r2unit_3*crsprod_2/mcrsprd;
    th2unit_2 = -crsprod_1/mcrsprd*r2unit_3 + r2unit_1*crsprod_3/mcrsprd;
    th2unit_3 = -crsprod_2/mcrsprd*r2unit_1 + r2unit_2*crsprod_1/mcrsprd;    % make 100.4% sure it's in (-1 <= x <= +1)
    
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
%     T0  = LancasterBlanchard(0, q, m);
    
    x_tmp = 0;
    q_tmp = q;
    m_tmp = m;
    
    LancasterBlanchard;
    T0 = T_tmp;
    
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
%             [dummy, Tp, Tpp, Tppp] = LancasterBlanchard(xM, q, m);%#ok
            x_tmp = xM;
            q_tmp = q;
            m_tmp = m;
            LancasterBlanchard;
            
            dummy = T_tmp;
            Tp = Tp_tmp;
            Tpp =Tpp_tmp;
            Tppp = Tppp_tmp;
            
            
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
        x_tmp = xM;
        q_tmp = q;
        m_tmp = m;
        
        
        LancasterBlanchard;
        TM = T_tmp;
        
%         TM = LancasterBlanchard(xM, q, m);
        
        % T should lie above the minimum T
        if (TM > T), exitflag = -1; return; end
        
        % find two initial values for second solution (again with lambda-type patch)
        % --------------------------------------------------------------------------
        
        % some initial values
        TmTM  = T - TM;   T0mTM = T0 - TM;
%         [dummy, Tp, Tpp] = LancasterBlanchard(xM, q, m);%#ok
        
        x_tmp = xM;
        q_tmp = q;
        m_tmp = m;
        
        LancasterBlanchard;
        dummy = T_tmp;
        Tp = Tp_tmp;
        Tpp = Tpp_tmp;
                
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
    x = x0; Tx = inf; iterations = 0;
    while abs(Tx) > tol        
        % iterations
        iterations = iterations + 1;        
        % compute function value, and first two derivatives
%         [Tx, Tp, Tpp] = LancasterBlanchard(x, q, m);
        x_tmp = x;
        q_tmp = q;
        m_tmp = m;
        
        LancasterBlanchard;
        Tx = T_tmp;
        Tp = Tp_tmp;
        Tpp = Tpp_tmp;
        
        % find the root of the *difference* between the
        % function value [T_x] and the required time [T]
        Tx = Tx - T;
        % new value of x
        xp = x;
        x  = x - 2*Tx*Tp ./ (2*Tp^2 - Tx*Tpp);                 
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
        sigma = 2*sqrt(r1*r2/(c^2)) * sin(dth/2)
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
    r1vec_1_tmp = r1vec_1;
    r1vec_2_tmp = r1vec_2;
    r1vec_3_tmp = r1vec_3;
    
    r1_tmp = r1;
    
%     r2vec_1_tmp = r2vec_1;
%     r2vec_2_tmp = r2vec_2;
%     r2vec_3_tmp = r2vec_3;

    r2vec_1_tmp = r1vec_1;
    r2vec_2_tmp = r1vec_2;
    r2vec_3_tmp = r1vec_3;
    
    r2_tmp = r2;
    
    dth_tmp = dth;
    a_tmp = a;
    V1_1_tmp = V1_1;
    V1_2_tmp = V1_2;
    V1_3_tmp = V1_3;
    V2_1_tmp = V2_1;
    V2_2_tmp = V2_2;
    V2_3_tmp = V2_3;
    
    m_tmp = m;
    muC_tmp = muC;
    
    minmax_distances;
    
    extremal_distances_L = extremal_distances_L_tmp;
    extremal_distances_H = extremal_distances_H_tmp;
    
%     extremal_distances = minmax_distances(r1vec, r1, r1vec, r2, dth, a, V1, V2, m, muC);
    
end