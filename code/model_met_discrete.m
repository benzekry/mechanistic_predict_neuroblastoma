function [x, nb, mass, mets_generation] = model_met_discrete(...
    T,...                   % final time
    dt,...                  % discretization step
    growth_model_PT,...     % PT growth rate - Volume of tumor at each timestep (including possible treatment)
    param_growth_PT,...     % PT growth parameters
    X0p,...                 % PT initial condition
    resection_time, ...     % time of removal of primary tumor
    growth_model_mets,...   % mets growth rate - volume
    param_growth_met,...    % mets growth parameters
    X0,...                  % mets initial condition (starting nb of cells to create a met)
    param_dissemination,... % dissemination parameters [mu, gamma, V_d0, tau, t_d0]
    flag_secondary_dissemination,...
    visible_s)              % object with fields related to simulation of only visible mets. flag = 1: simulate only visible mets. If 1, then need fields and visible_threshold (diameter in mm)
%--------------------------------------------------------------------------
% Implementation of a model for metastatic growth (dissemination and colonization)
% This code considers unperturbed growth only (autonomous case, i.e no therapy)
% Only one characteristic is generated at the beginning and reused all
% along the algorithm
% Discrete (non-stochastic) emission of metastases (expectation of a Poisson process)
% Includes the possibility of a delay to metastatic initiation, dormancy
% and secondary dissemination (metastases from metastases)
%--------------------------------------------------------------------------
mu    = param_dissemination(1);
gamma = param_dissemination(2);
if length(param_dissemination) >= 4
	V_d0  = param_dissemination(3);
	tau0  = param_dissemination(4);    
else
	V_d0  = 0;
	tau0  = 0;    
end
if length(param_dissemination) >= 5
    t_d0  = param_dissemination(5);
else
    t_d0  = 0;
end
M     = length(X0);
%--------------------------------------------------------------------------
% Discretization parameters
%--------------------------------------------------------------------------
temps = 0:dt:T;
K     = length(temps)-1;
%--------------------------------------------------------------------------
% Primary tumor growth
%--------------------------------------------------------------------------
Xp = growth_model_PT(param_growth_PT, temps, X0p); % PT volume with treatment
Xp = reshape(Xp, K+1, 1);
    function y = emission_rate(V)
        if V > 0
            y = mu*V.^gamma;
        else 
            y = 0;
        end
    end
%--------------------------------------------------------------------------
% Builds L(t) = \int_0^t d(V_p(s))ds so that it is not recomputed each
% time
%--------------------------------------------------------------------------
L = zeros(1, K+1); % an array with the nb of mets at each time step
for k = 2:K+1
    tk_moins_un = temps(k-1);
    tk          = temps(k);
    deltat      = tk-tk_moins_un;
    L(k)        = L(k-1) + deltat/2*(emission_rate(Xp(k-1))+emission_rate(Xp(k)));
end
%--------------------------------------------------------------------------
% Initialization
%--------------------------------------------------------------------------
x                = cell(1, K+1);  % 1-by-K+1 array of empty matrices
nb               = zeros(K+1, M);
mass             = zeros(K+1, M);
emission_indices = [];
accumulated      = 0;
accumulated_mets = [];
mets_generation  = [];
if ~exist('flag_secondary_dissemination', 'var')
   flag_secondary_dissemination = 0;
end
%options = optimset('TolX', 0.01);
%--------------------------------------------------------------------------
% Characteristic (only one)
%--------------------------------------------------------------------------
if tau0 > 0
    X                    = zeros(size(temps));
    idx_dorm             = find((temps - tau0) <= 0);
    idx_no_dorm          = find((temps - tau0) > 0);
    X(idx_dorm)          = X0*ones(size(temps(idx_dorm)));
    temps_shift          = temps - tau0;
    X(idx_no_dorm)       = growth_model_mets(param_growth_met, temps_shift(idx_no_dorm), X0);  % growth of 1 met from the size X0
else
    X = growth_model_mets(param_growth_met, temps, X0);  % growth of 1 met from the size X0
end
X = reshape(X, K+1, 1);
%--------------------------------------------------------------------------
% Time to visible size
%--------------------------------------------------------------------------
function tau = time_to_size_met(S)
    idx_tau = find(X > S, 1, 'first');
    tau     = temps(idx_tau);
end
if visible_s.flag == 1
    visible_threshold_cell = diam2cell(visible_s.visible_threshold);    
    tau_vis                = time_to_size_met(visible_threshold_cell);
    if isempty(tau_vis) == 1
        tau_vis = inf;
    end
end
%--------------------------------------------------------------------------
% Algorithm
%--------------------------------------------------------------------------
Xk_minus_one = X(2);
t_d0_idx     = max(find(Xp >= V_d0, 1, 'first'), find(temps > t_d0, 1, 'first')); 
for k = t_d0_idx+1:K+1    
    tk = temps(k);
    if mod(k, 400) == 0
        prc  =  k/K*100;
%         display(['avancement  =  ' num2str(prc)]);        
    end
    nb_new = 0;    
    %----------------------------------------------------------------------
    % Mets growth
    %----------------------------------------------------------------------
    Xk               = X(k - emission_indices + 1);  % size of mets at time tk
    mets_nb          = length(Xk);
    %----------------------------------------------------------------------
    % Emission of new mets during dt (=\int_t^{t+dt} d(V_p(s))ds)   
    % * If restriction to only visible, then emission only if time is smaller
    % T - tau_vis (if not, it won't be visible at time T)
    % * Only if time is smaller than resection time
    %----------------------------------------------------------------------    
    if (visible_s.flag == 0) || (tk < (T - tau_vis))        
        if tk < resection_time
            accumulated      = accumulated + ...
                dt/2*(emission_rate(Xp(k-1)) + emission_rate(Xp(k)));
            if accumulated >= 1
                nb_new                = round(accumulated); % new met(s)
                accumulated           = accumulated - nb_new;
                Xk                    = [Xk; X0.*ones(nb_new, M)];
                for count = 1:nb_new   % if more than one met have been established
                    accumulated_mets      = [accumulated_mets; 0];
                    mets_generation       = [mets_generation;  1]; % mets from PT are first generation
                end
            else
                nb_new                = 0;
            end
        end
        %----------------------------------------------------------------------
        % Emission of secondary mets
        %----------------------------------------------------------------------
        if flag_secondary_dissemination == 1
            %------------------------------------------------------------------
            % Calculation of emission by each existing met
            %------------------------------------------------------------------
            for j = 1:mets_nb
                nb_new_j = 0;
                % to add a delay effect add here an a condition
                if (Xk(j) >= V_d0)
                    if length(Xk_minus_one) < j
                      Xk_minus_one = [Xk_minus_one; 0];
                    end            
                    accumulated_mets(j) = accumulated_mets(j) + ...
                    dt/2*(emission_rate(Xk(j)) + emission_rate(Xk_minus_one(j)));            
                    if accumulated_mets(j) >= 1
                       nb_new_j            = round(accumulated_mets(j)); % new met(s)
                       accumulated_mets(j) = accumulated_mets(j) - nb_new_j;
                       Xk                  = [Xk; X0.*ones(nb_new_j, M)];
                       for count = 1:nb_new_j         % if more than one met have been established
                          accumulated_mets = [accumulated_mets; 0];
                          mets_generation  = [mets_generation; mets_generation(j) + 1];
                       end                   
                    end 
                end
                nb_new = nb_new + nb_new_j;
            end
        end    
    end
    emission_indices     = [emission_indices, k*ones(1,nb_new)];
    x{k}                 = Xk;
    nb(k)                = length(Xk);
    mass(k)              = sum(Xk);    
    Xk_minus_one         = Xk;
end
function N = diam2cell(D)
% Converts diameter in mm into number of cells
N = 1/6*pi*D.^3*1e6;
end
end