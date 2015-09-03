clear
clc
mu = 1.327124e11; % Sun Gravitational parameter (km^3/s^2) 

planet_id_i = 3; % Planet identification 3---> Earth
planet_id_f = 4; % Planet identification 4---> Mars

% Initial departure date
year = 2005;
month = 1;
day = 1 ;

% Initial Julian date
JD_0 = J0(year, month ,day);

% Create a vector for the departure Julian date over a compelete year with
% increments of one day
JD_d_vec = JD_0:1:JD_0 + 1*365;
Length_JD_d_vec = length(JD_d_vec);

% Create a vector for the time of flight with increments of 5 days
TOF_vec = 100:5:2*365;
Length_TOF_vec = length(TOF_vec);

% Allocate variable to store the position and velocity vectors for each of
% the Julian dates and time of flight
GridDataR_i = zeros(Length_JD_d_vec,3);
GridDataV_i = zeros(Length_JD_d_vec,3);
GridDataR_f = zeros(Length_TOF_vec,3);
GridDataV_f = zeros(Length_TOF_vec,3);

GPU_Ri_1_tmp = [];
GPU_Ri_2_tmp = [];
GPU_Ri_3_tmp = [];

GPU_Rf_1_tmp = [];
GPU_Rf_2_tmp = [];
GPU_Rf_3_tmp = [];

GPU_TOF_tmp = [];

for Jd_index=1:Length_JD_d_vec; % loop over the departure Julian date
    JD_d = JD_d_vec(Jd_index);
    [coe_i, R_i, V_i, jd_i] = planet_elements_and_sv(planet_id_i,JD_d, mu); % calculate the position and velocity  
    GridDataR_i(Jd_index,:) = R_i;
    GridDataV_i(Jd_index,:) = V_i;
    for TOF_index = 1:Length_TOF_vec % loop over the time of fliht at each departure date
        TOF = TOF_vec(TOF_index);
        JD_f = JD_d + TOF;
        [coe_f, R_f, V_f, jd_f] = planet_elements_and_sv(planet_id_f,JD_f, mu); %calculate the position and velocity  
        GridDataR_f(TOF_index,:) = R_f;
        GridDataV_f(TOF_index,:) = V_f;
        
        GPU_Ri_1_tmp = [GPU_Ri_1_tmp; R_i(1)];
        GPU_Ri_2_tmp = [GPU_Ri_2_tmp; R_i(2)];
        GPU_Ri_3_tmp = [GPU_Ri_3_tmp; R_i(3)];
        
        GPU_Rf_1_tmp = [GPU_Rf_1_tmp; R_f(1)];
        GPU_Rf_2_tmp = [GPU_Rf_2_tmp; R_f(2)];
        GPU_Rf_3_tmp = [GPU_Rf_3_tmp; R_f(3)];
        
        GPU_TOF_tmp = [GPU_TOF_tmp; TOF];
        
        
    end
end


GPU_Ri_1 = gpuArray(GPU_Ri_1_tmp);
GPU_Ri_2 = gpuArray(GPU_Ri_2_tmp);
GPU_Ri_3 = gpuArray(GPU_Ri_3_tmp);

GPU_Rf_1 = gpuArray(GPU_Rf_1_tmp);
GPU_Rf_2 = gpuArray(GPU_Rf_2_tmp);
GPU_Rf_3 = gpuArray(GPU_Rf_3_tmp);

GPU_TOF = gpuArray(GPU_TOF_tmp);

% GridDataR_f_GPU = gpuArray(GridDataR_f);
% 
% GridDataR_i_1 = GridDataR_i(1:127,1);
% GridDataR_i_2 = GridDataR_i(1:127,2);
% GridDataR_i_3 = GridDataR_i(1:127,3);
% 
% GridDataR_f_1 = GridDataR_f(1:127,1);
% GridDataR_f_2 = GridDataR_f(1:127,2);
% GridDataR_f_3 = GridDataR_f(1:127,3);

num_orbit_vec = ones([length(GPU_Ri_1),1])*0;
mu_vec = ones([length(GPU_Ri_1),1])*mu;

tic
[V1_1, V1_2, V1_3, V2_1, V2_2, V2_3, extremal_distances_L, extremal_distances_H, exitflag] ...
    = arrayfun(@lambert_GPU, GPU_Ri_1, GPU_Ri_2, GPU_Ri_3, GPU_Rf_1, GPU_Rf_2, GPU_Rf_3, GPU_TOF, num_orbit_vec, mu_vec);
toc

V1_final = [gather(V1_1) gather(V1_2) gather(V1_3)];
V2_final = [gather(V2_1) gather(V2_2) gather(V2_3)];
extremal_distances_final = [gather(extremal_distances_L) gather(extremal_distances_H)];
exitflag_final = gather(exitflag);


