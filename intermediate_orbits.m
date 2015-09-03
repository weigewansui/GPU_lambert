function [rECI,vECI,t_nodes] = intermediate_orbits(coe_i,coe_f,N_rev,N_nodes,tf)
%--------------------------------------------------------------------------
ei    = coe_i(2);
RAi   = coe_i(3);
incli = coe_i(4);
wi    = coe_i(5);
TAi   = coe_i(6);
ai    = coe_i(7);
Pi    = ai*(1-ei^2);
%--------------------------------------------------------------------------
ef    = coe_f(2);
RAf   = coe_f(3);
inclf = coe_f(4);
wf    = coe_f(5);
TA2   = coe_f(6);
af    = coe_f(7);
pf    = af*(1-ef^2);
%--------------------------------------------------------------------------
TAf = TA2+2*pi*(N_rev+1);
[ae,be] = linear_coe_coeffs(TAi,TAf,ei,ef);
[aW,bW] = linear_coe_coeffs(TAi,TAf,RAi,RAf);
[ai,bi] = linear_coe_coeffs(TAi,TAf,incli,inclf);
[aw,bw] = linear_coe_coeffs(TAi,TAf,wi,wf);
[ap,bp] = linear_coe_coeffs(TAi,TAf,Pi,pf);
[at,bt] = linear_coe_coeffs(0  ,tf ,TAi,TAf);
%--------------------------------------------------------------------------
t_nodes = linspace(0,tf,N_nodes);
theta = at*t_nodes+bt;
e = ae*theta+be;
W = aW*theta+bW;
incl = ai*theta+bi;
w = aw*theta+bw;
p = ap*theta+bp;
%--------------------------------------------------------------------------
length_t = length(theta);
rECI = zeros(3,length_t);
vECI = zeros(3,length_t);
for i=1:length_t
    r = p(i)/(1+e(i)*cos(theta(i)));
    rp = [r*cos(theta(i));r*sin(theta(i));0];
    vp = [-sin(theta(i));e(i)+cos(theta(i));0]/sqrt(p(i));
    R3_W = [ cos(W(i)) sin(W(i)) 0
        -sin(W(i)) cos(W(i)) 0
        0       0  1];
    %...Equation 4.40:
    R1_i = [1         0         0
        0  cos(incl(i)) sin(incl(i))
        0 -sin(incl(i)) cos(incl(i))];
    %...Equation 4.41:
    R3_w = [ cos(w(i)) sin(w(i)) 0
            -sin(w(i)) cos(w(i)) 0
        0 0 1];
    TP2ECI = R3_W.'*R1_i.'*R3_w.';
    rECI(:,i) = TP2ECI*rp;
    vECI(:,i) = TP2ECI*vp;
end

% plot(t_nodes,rECI(1,:));
end