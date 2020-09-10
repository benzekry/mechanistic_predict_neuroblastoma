function birth_times = simulate_f(idx_iter)
table_res       = readtable('mechanistic/df_fit_indiv.csv');
patient_id      = table_res{idx_iter, 'id'};
folder_glob = ['simulations/patient_' num2str(patient_id)];
if ~exist(folder_glob, 'dir')
    mkdir(folder_glob)
end
%-----------------------------------------------
%% Parameters
%------------------------------------------------------------------------
median_DT       = 48; % (hours)
% Model parameters
visible_threshold_SIOPEN = table_res{idx_iter, 'visible_threshold'};
V_diag          = vol2cell(table_res{idx_iter, 'tumor_size_mm3'}); % (cell)
alpha           = log(2)/median_DT*24; % (day^-1)
mu              = table_res{idx_iter, 'mu'}; % (cell^-1.day^-1)
% More generic
Vc              = 1; % (cell)
V0              = 1; % (cell)
dt              = 0.01; % (day)
Tmax            = 365; % (day)
temps_Vp        = 0:dt:Tmax; % (day)
growth_model    = @(param, time, X0) X0*exp(param*time);
Vp_loc          = growth_model(alpha, temps_Vp, V0);
idx_diag        = find(Vp_loc > V_diag, 1, 'first');
T_diag          = temps_Vp(idx_diag);
resection_time  = T_diag;
T_end           = T_diag;
secondary_diss  = 0;
detection_limit = diam2cell(3);
temps           = 0:dt:(T_end + dt); % total simulation time
% data_s
data_s                = struct();
data_s.T1             = T_diag;
data_s.Ts             = T_end;
data_s.vis_thresh     = 3; % (mm)
data_s.dataset_name   = 'children_neuroblastoma';
data_s.resection_flag = 0;
data_s.T_resec_date   = []; % date of resection of mets
% model_s
model_s                         = struct();
model_s.growth_model_PT         = growth_model;
model_s.growth_model_met        = growth_model;
model_s.S0p                     = 1;
model_s.S0                      = 1;
model_s.param_s.names           = {'alpha', 'mu'};
model_s.param_s.param_all       = [alpha, mu];
model_s.param_s.mets_growth_idx = [1];
model_s.param_s.PT_all          = alpha;
model_s.unit_met_data2model     = @diam2cell;
model_s.unit_met_model2data     = @cell2diam;
model_s.param_s.PT_ther         = [];
% visible_s. to define if computing only visible or not
visible_s.flag                  = 1;
visible_s.visible_threshold     = visible_threshold_SIOPEN;
%------------------------------------------------------------------------
%% Cosmetics
%------------------------------------------------------------------------
flag_logX                        = 0;
hist_movie                       = 0;
% flag_s
flag_s.add_data_to_all_log_kinetics = 0;
flag_s.secondary_dissemination   = 0;
% output_s
output_s                         = struct();
output_s.post_diag               = 1;
output_s.plot_mets_distribution  = 0;
output_s.plot_mets_distrib_final = 0;
output_s.time_unit               = 30;
output_s.time_unit_name          = 'months';
output_s.time_label_nb_vis       = 'Time (months)';
output_s.size_range              = 1:10:40;
output_s.xmax                    = 40;
output_s.ymax                    = 20;
output_s.unit                    = 'a.u.'; % unit label
output_s.plot_only_visible       = 0;
output_s.plot_birth_times        = 0;
output_s.manips_birth_times.ymax_mets = 100;
output_s.manips_birth_times.vis_thresh_SIOPEN = visible_threshold_SIOPEN;
%------------------------------------------------------------------------
%% Plot PT
%------------------------------------------------------------------------
% Vp              = growth_model(model_s.param_s.PT_all, temps, V0);
% semilogy(temps/output_s.time_unit, Vp)
%------------------------------------------------------------------------
%% Simulate mets
%------------------------------------------------------------------------
[x, ~, ~, generation]  = model_met_discrete(...
    T_end+dt,...
    dt,...
    growth_model,...
    [alpha],...
    V0,...
    T_diag,...
    growth_model,...
    [alpha],...
    V0,...
    [mu, 1],...
    secondary_diss,...
    visible_s);
save([folder_glob '/simu_mets'], 'x')
save([folder_glob '/fit']);
%------------------------------------------------------------------------
%% Plot mets distributions
%------------------------------------------------------------------------
Vp    = growth_model(model_s.param_s.PT_all, temps, V0);
Vp    = Vp.*(temps < T_diag);
if output_s.plot_mets_distribution == 1
    plot_mets_distribution(...
        temps,...
        x,...
        idx_diag,...
        folder_glob,...
        output_s.unit,...
        Vp,...
        detection_limit,... % (cell)
        flag_logX,...
        output_s.xmax,...
        output_s.ymax,...
        hist_movie,...
        data_s,...
        output_s)
    close all
end
%------------------------------------------------------------------------
%% Plot birth times
%------------------------------------------------------------------------
birth_times = manips_birth_times(...
    folder_glob,...
    temps,...
    x,...
    data_s,...
    model_s,...
    flag_s,...
    output_s);
close all
function N = vol2cell(V)
% Converts volume in mm3 into number of cells
N = V*1e6;
function N = diam2cell(D)
% Converts diameter in mm into number of cells
N = 1/6*pi*D.^3*1e6;
