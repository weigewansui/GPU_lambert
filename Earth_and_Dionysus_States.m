function [RE,VE,RD,VD]= Earth_and_Dionysus_States(~)
clc
AU=1.49597870691e8;
muS=1.32712440018e11;
dtr=pi/180;
% t=2459832.519607843400-2400000.5; % MJD at required time
t = 51543+4741+3534; %51543 is for Jan 2000 and 4741 is Departure date (MJD 2000) and 3534 is Time of flight 
%-------------------------- Dionysus Orbit Date ------------------------------
t0E = 53400; %Epoch Date
aE = 2.2*AU;
eE = 0.542;
inclE = 13.6*dtr;
OmegaE = 82.2*dtr;
wE = 204.2*dtr;
M0E = 114.4232*dtr;
ME = M0E+sqrt(muS/aE^3)*(t-t0E)*86400;
ME = wrapToPi(ME);
EE = kepler_Eq(eE,ME);
thetaE = 2*atan(sqrt((1+eE)/(1-eE))*tan(EE/2));
rE = aE*(1-eE^2)/(1+eE*cos(thetaE));
gamaE = atan((eE*sin(thetaE))/(1+eE*cos(thetaE)));
vE = sqrt(2*muS/rE-muS/aE);
%------------------------ Earth position and velocity ---------------------
xD = rE*( cos(thetaE+wE)*cos(OmegaE)-sin(thetaE+wE)*cos(inclE)*sin(OmegaE));
yD = rE*( cos(thetaE+wE)*sin(OmegaE)+sin(thetaE+wE)*cos(inclE)*cos(OmegaE));
zD = rE*( sin(thetaE+wE)*sin(inclE));
RD = [xD;yD;zD];
%--------------------------------------------------------------------------
vxD= vE*(-sin(thetaE+wE-gamaE)*cos(OmegaE)-cos(thetaE+wE-gamaE)*cos(inclE)*sin(OmegaE));
vyD= vE*(-sin(thetaE+wE-gamaE)*sin(OmegaE)+cos(thetaE+wE-gamaE)*cos(inclE)*cos(OmegaE));
vzD= vE*( cos(thetaE+wE-gamaE)*sin(inclE));
VD = [vxD;vyD;vzD];
%-------------------------- Earth Orbit Date ------------------------------
t = 51543+4741;
t0E=54000;
aE=0.999988049532578*AU;
eE=1.671681163160e-2;
inclE=0.8854353079654e-3*dtr;
OmegaE=175.40647696473*dtr;
wE=287.61577546182*dtr;
M0E=257.60683707535*dtr;
ME=M0E+sqrt(muS/aE^3)*(t-t0E)*86400;
ME=wrapToPi(ME);
EE=kepler_Eq(eE,ME);
thetaE=2*atan(sqrt((1+eE)/(1-eE))*tan(EE/2));
rE=aE*(1-eE^2)/(1+eE*cos(thetaE));
gamaE=atan((eE*sin(thetaE))/(1+eE*cos(thetaE)));
vE=sqrt(2*muS/rE-muS/aE);
%------------------------ Earth position and velocity ---------------------
xE = rE*( cos(thetaE+wE)*cos(OmegaE)-sin(thetaE+wE)*cos(inclE)*sin(OmegaE));
yE = rE*( cos(thetaE+wE)*sin(OmegaE)+sin(thetaE+wE)*cos(inclE)*cos(OmegaE));
zE = rE*( sin(thetaE+wE)*sin(inclE));
RE = [xE;yE;zE];
%--------------------------------------------------------------------------
vxE= vE*(-sin(thetaE+wE-gamaE)*cos(OmegaE)-cos(thetaE+wE-gamaE)*cos(inclE)*sin(OmegaE));
vyE= vE*(-sin(thetaE+wE-gamaE)*sin(OmegaE)+cos(thetaE+wE-gamaE)*cos(inclE)*cos(OmegaE));
vzE= vE*( cos(thetaE+wE-gamaE)*sin(inclE));
VE = [vxE;vyE;vzE];
end