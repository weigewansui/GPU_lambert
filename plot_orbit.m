function [rECI,rpe_ECI,rap_ECI] = plot_orbit(coeE)
% mus = 1.32712440018e+11; %Km^3/sec^2
mus = 1;
h = coeE(1);
e = coeE(2);
RA = coeE(3);
incl = coeE(4);
w = coeE(5);
TA = coeE(6);
a = coeE(7);
theta = linspace(TA,TA+2*pi,200);
r = h^2/mus*1./(1+e*cos(theta));
rp = [r.*cos(theta);r.*sin(theta);zeros(1,length(theta))];
R3_W = [ cos(RA) sin(RA) 0
        -sin(RA) cos(RA) 0
              0       0  1];
%...Equation 4.40:
R1_i = [1         0         0
        0  cos(incl) sin(incl)
        0 -sin(incl) cos(incl)];
%...Equation 4.41:
R3_w = [ cos(w) sin(w) 0
        -sin(w) cos(w) 0
            0 0 1];
TP2ECI = R3_W'*R1_i'*R3_w';
rECI = TP2ECI*rp;
rpe = [a*(1-e);0;0];
rpe_ECI = TP2ECI*rpe;
rap = [-a*(1+e);0;0];
rap_ECI = TP2ECI*rap;
end