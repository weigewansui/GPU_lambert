% function [T, Tp, Tpp, Tppp] = LancasterBlanchard(x, q, m)
    
    % protection against idiotic input
    if (x_tmp < -1) % impossible; negative eccentricity
        x_tmp = abs(x_tmp) - 2;
    elseif (x_tmp == -1) % impossible; offset x slightly
        x_tmp = x_tmp + eps;
    end
    
    % compute parameter E
    E  = x_tmp*x_tmp - 1;    
    
    % T(x), T'(x), T''(x)
    if x_tmp == 1 % exactly parabolic; solutions known exactly
        % T(x)
        T_tmp = 4/3*(1-q_tmp^3);
        % T'(x)
        Tp_tmp = 4/5*(q_tmp^5 - 1);
        % T''(x)
        Tpp_tmp = Tp_tmp + 120/70*(1 - q_tmp^7);
        % T'''(x)
        Tppp_tmp = 3*(Tpp_tmp - Tp_tmp) + 2400/1080*(q_tmp^9 - 1);
        
    elseif abs(x_tmp-1) < 1e-2 % near-parabolic; compute with series
        % evaluate sigma
%         [sig1, dsigdx1, d2sigdx21, d3sigdx31] = sigmax(-E);
        y_tmp = -E;
        sigm;
        sig1 = sig;
        dsigdx1 = dsigdx;
        d2sigdx21 = d2sigdx2;
        d3sigdx31 = d3sigdx3;
        
%         [sig2, dsigdx2, d2sigdx22, d3sigdx32] = sigmax(-E*q*q);  
        
        y_tmp = -E*q_tmp*q_tmp;
        sigm;
        sig2 = sig;
        dsigdx2 = dsigdx;
        d2sigdx22 = d2sigdx2;
        d3sigdx32 = d3sigdx3;
        
        % T(x)
        T_tmp = sig1 - q_tmp^3*sig2;
        % T'(x)
        Tp_tmp = 2*x_tmp*(q_tmp^5*dsigdx2 - dsigdx1);
        % T''(x)        
        Tpp_tmp = Tp_tmp/x_tmp + 4*x_tmp^2*(d2sigdx21 - q_tmp^7*d2sigdx22);
        % T'''(x)
        Tppp_tmp = 3*(Tpp_tmp-Tp_tmp/x_tmp)/x_tmp + ...
            8*x_tmp*x_tmp*(q_tmp^9*d3sigdx32 - d3sigdx31);
        
    else % all other cases
        % compute all substitution functions
        y_inside  = sqrt(abs(E));
        z  = sqrt(1 + q_tmp^2*E);
        f  = y_inside*(z - q_tmp*x_tmp);
        g  = x_tmp*z - q_tmp*E;

        % BUGFIX: (Simon Tardivel) this line is incorrect for E==0 and f+g==0
        % d  = (E < 0)*(atan2(f, g) + pi*m) + (E > 0)*log( max(0, f + g) );
        % it should be written out like so: 
        if (E<0)
            d = atan2(f, g) + pi*m_tmp;
        elseif (E==0)
            d = 0;
        else 
            d = log(max(0, f+g));
        end

        % T(x)
        T_tmp = 2*(x_tmp - q_tmp*z - d/y_inside)/E;
        %  T'(x)
        
        Tp_tmp = (4 - 4*q_tmp^3*x_tmp/z - 3*x_tmp*T_tmp)/E;
        % T''(x)
        Tpp_tmp = (-4*q_tmp^3/z * (1 - q_tmp^2*x_tmp^2/z^2) - 3*T_tmp - 3*x_tmp*Tp_tmp)/E;
        % T'''(x) 
        Tppp_tmp = (4*q_tmp^3/z^2*((1 - q_tmp^2*x_tmp^2/z^2) + 2*q_tmp^2*x_tmp/z^2*(z - x_tmp)) - 8*Tp_tmp - 7*x_tmp*Tpp_tmp)/E  ; 
        
    end
% end