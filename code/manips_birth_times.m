function birth_times = manips_birth_times(...
    folder,...
    temps,...
    x,...
    data_s,...
    model_s,...
    flag_s,...
    output_s)
%--------------------------------------------------------------------------------------------------
% Plots birth times and growth curves from fit results of the metastatic model to size distributions
% data
%--------------------------------------------------------------------------------------------------
if ~exist('data_s', 'var')
    res      = load([folder '/fit.mat']);
    temps    = res.temps;
    x        = res.x;
    data_s   = res.data_s;
    model_s  = res.model_s;
    flag_s   = res.flag_s;
    output_s = res.output_s;
end
if ~isfield(output_s, 'plot_birth_times')
    output_s.plot_birth_times = 1
end
if isfield(model_s, 'fit')
    param_dict  = containers.Map(model_s.param_s.names, model_s.fit.param_all);
    param_all   = model_s.fit.param_all;
else
    param_dict  = containers.Map(model_s.param_s.names, model_s.param_s.param_all);
    param_all   = model_s.param_s.param_all;
end
T           = data_s.Ts(end);
idx         = find(temps > T, 1, 'first');
vis_thresh  = model_s.unit_met_data2model(data_s.vis_thresh);
x_vis_loc   = x{idx}(x{idx} > vis_thresh);
N_vis       = length(x_vis_loc);
N           = length(x{end});
if ~isfield(output_s, 'plot_only_visible')
    N_plot = N_vis;
elseif output_s.plot_only_visible == 0
    N_plot = N;
end
birth_times = zeros(1, N);
K           = length(temps);
comp        = 1;
scale       = output_s.time_unit;
colors      = get(gca, 'ColorOrder');
%--------------------------------------------------------------------------
% Birth times
%--------------------------------------------------------------------------
for k = 1:K
    if length(x{k}) >= comp
        birth_times(comp) = temps(k);
        comp              = comp +1;
    end
end
% Plot birth times
if output_s.post_diag == 1
    shift_x = data_s.T1;
else
    shift_x = 0;
end
if output_s.plot_birth_times == 1
    t_max = (T - shift_x)/scale + 1;
    figure(1)
    clf
    plot([-shift_x/scale, (data_s.Ts(end)-shift_x)/scale], [0, 0])
    set(gca, 'XLim', [-shift_x/scale, (data_s.Ts(end)-shift_x)/scale])
    set(gca, 'YLim', [-0.5, 0.5])
    hold on
    for n = 1:N_plot
        line(...
            [(birth_times(n) - shift_x)/scale (birth_times(n) - shift_x)/scale],...
            [-0.1, 0.1],...
            'Color', 'k')
    end
    hold off
    %--------------------------------------------------------------------------
    % Line at diagnosis and delay, if any
    %--------------------------------------------------------------------------
    y_lim = get(gca, 'Ylim');
    line([(data_s.T1 - shift_x)/scale, (data_s.T1 - shift_x)/scale], y_lim, 'Linestyle', '--', 'Color', 'k')
    if isKey(param_dict, '$t_d$') && param_dict('$t_d$') > 0
        line([(param_dict('$t_d$') - shift_x)/scale, (param_dict('$t_d$') - shift_x)/scale], y_lim, 'Linestyle', '--', 'Color', 'r')
    end
    %--------------------------------------------------------------------------
    % Labels, cosmetics and export
    %--------------------------------------------------------------------------
    xlabel(output_s.time_label_nb_vis)
    set_fonts_lines(gca)
    set(gca, 'ytick', [])
    set(gca, 'yticklabel', [])
    export_fig([folder '/birth_times.pdf']);
end
%--------------------------------------------------------------------------
% Growth curves with primary tumor
%--------------------------------------------------------------------------
figure(2)
clf
% Generate two axes
[AX, H1, H2] = plotyy(1, 1, 1, 1);
set(H1, 'LineStyle', 'none', 'Marker', 'none');
set(H2, 'LineStyle', 'none', 'Marker', 'none');
xlim = [-shift_x/scale, (data_s.Ts(end)-shift_x)/scale];
set(AX, {'XLim'}, {xlim; xlim})
set(AX, {'Yscale'}, {'log'; 'log'})
set(AX(1), 'xlim', [-3, 0.5])
set(AX(2), 'xlim', [-3, 0.5]);
set(AX(1), 'ylim', [1, 1e13])
set(AX(2), 'ylim', [cell2diam(1), cell2diam(1e13)]);
set(AX, {'ycolor'}, {'k';'k'});
set(AX, {'Xticklabel'}, {[];[]});
set(AX, {'Xtick'}, {[];[]});
set(AX, {'Yticklabel'}, {[];[]});
set(AX, {'Ytick'}, {[];[]});
hold(AX(1), 'on')
hold(AX(2), 'on')
axes(AX(1))
% Plot PT
Vp = model_s.growth_model_PT(...
    model_s.param_s.PT_all,...
    temps,...
    model_s.S0p);
% set to zero after surgery, if any
if isfield(data_s, 'resection_flag') && (data_s.resection_flag == 1)
    Vp = Vp.*(temps <= data_s.T1);
end
plot((temps - shift_x)/scale, Vp, 'color', colors(1, :))
hold on
if flag_s.add_data_to_all_log_kinetics == 1
    % Data
    plot(data_s.time_PT{1}/scale, vol2cell(data_s.PT{1}*1e3), 'o', 'color', colors(1, :))
end
if strfind(data_s.dataset_name, 'mice')
    hold off
end
% Metastases
for n = 1:N_plot
    if (flag_s.add_data_to_all_log_kinetics == 1) && (n <= length(data_s.mets))
        % Data
        plot(AX(1), data_s.time_mets{n}, diam2cell(data_s.mets{n}), 'or')
    end
    % Model
    if isKey(param_dict, '$\tau$') && (param_dict('$\tau$') > 0)
        %------------------------------------------------------------------
        % Dormancy
        %------------------------------------------------------------------
        tau       = param_dict('$\tau$');
        temps_loc = temps(temps >= (birth_times(n) + tau));
        V         = model_s.growth_model_met(...
            param_all(model_s.param_s.mets_growth_idx),...
            temps_loc - (birth_times(n) + tau),...
            model_s.S0);
        plot(AX(1), [(birth_times(n) - shift_x)/scale, (birth_times(n) + tau - shift_x)/scale], [1, 1], 'r')
    else
        %------------------------------------------------------------------
        % No dormancy
        %------------------------------------------------------------------
        temps_loc = temps(temps >= birth_times(n));
        V         = model_s.growth_model_met(...
            param_all(model_s.param_s.mets_growth_idx),...
            temps_loc - birth_times(n),...
            model_s.S0);
    end
    plot(AX(1), (temps_loc - shift_x)/scale, V, 'r')
end
hold off
set(AX(1), 'Yscale', 'log')
% if isfield(output_s, 'y_scale_size')
%     set(gca, 'Yscale', output_s.y_scale_size)
% else
%     set(gca, 'Yscale', 'log')
%     set(gca, 'ylim', [1 1e12])
% end
% set(gca, 'XLim', [-shift_x/scale, (data_s.Ts(end)-shift_x)/scale])
%--------------------------------------------------------------------------
% Vertical lines at diagnosis and delay, if any
%--------------------------------------------------------------------------
y_lim = get(gca(), 'Ylim');
hold(AX(1), 'on')
plot(AX(1), [(data_s.T1 - shift_x)/scale, (data_s.T1 - shift_x)/scale], y_lim, 'Linestyle', '--', 'Color', 'k')
if isKey(param_dict, '$t_d$') && param_dict('$t_d$') > 0
    line([(param_dict('$t_d$') - shift_x)/scale, (param_dict('$t_d$') - shift_x)/scale], y_lim, 'Linestyle', '--', 'Color', 'r')
end
hold off
%--------------------------------------------------------------------------
% Vertical line at TTR, if any
%--------------------------------------------------------------------------
if isfield(data_s, 'TTR')
    line(AX(1), [(data_s.T1 + data_s.TTR - shift_x)/scale, (data_s.T1 + data_s.TTR - shift_x)/scale],...
        y_lim,...
        'Linestyle', '--',...
        'Color', 'r')
end
%--------------------------------------------------------------------------
% Horizontal line at visible threshold
%--------------------------------------------------------------------------
vis_thresh_loc = diam2cell(output_s.manips_birth_times.vis_thresh_SIOPEN);
hold on
x_lim = get(gca, 'xlim');
text(x_lim(1), vis_thresh_loc*8, 'Detection limit')
line(get(gca, 'xlim'), [vis_thresh_loc, vis_thresh_loc], 'Linestyle', '--', 'Color', 'k')
%--------------------------------------------------------------------------
% Labels, cosmetics and export
%--------------------------------------------------------------------------
fontaxes = 20;
hold off
xlabel(output_s.time_label_nb_vis)
ylabel(['Tumor size (cells)'])
hold(AX(1), 'off')
hold(AX(2), 'off')
set(AX(1), 'XtickMode', 'auto')
set(AX(1), 'XticklabelMode', 'auto')
set(AX(1), 'YtickMode', 'auto')
set(AX(1), 'YtickLabelMode', 'auto')
set(AX(1), 'box', 'off')
set(AX(2), 'YtickMode', 'auto');
set(AX(2), 'YtickLabelMode', 'auto')
set(AX(1), 'Fontsize', fontaxes)
set(AX(2), 'Fontsize', fontaxes)
xlabel(AX(1), output_s.time_label_nb_vis)
ylabel(AX(1), ['Tumor size (cells)'])
ylabel(AX(2), ['Tumor size (mm)'])
AX(1).Position(3) = AX(1).Position(3)*0.95; % to prevent axis labels to be cut
AX(1).Position(2) = AX(1).Position(2)*1.1; % to prevent axis labels to be cut
set_fonts_lines(AX(1))
set_fonts_lines(AX(2))
AX(2).YTickLabel = AX(2).YTick;
export_fig([folder '/growth_mets.pdf']);
%--------------------------------------------------------------------------
% Growth curves only metastases
%--------------------------------------------------------------------------
if isfield(data_s, 'met')
    folder_kinetics = [folder '/mets_kinetics'];
    if ~exist(folder_kinetics, 'dir')
        mkdir(folder_kinetics)
    end
    figure(3)
    hold on
    for idx_met = 1:length(data_s.mets)
        % Data
        met = data_s.mets{idx_met};
        figure(3)
        errorbar(data_s.time_mets{idx_met}, met, 0.1*met, 'o', 'Color', colors(1, :))
        figure(4) % plot of each met
        errorbar(data_s.time_mets{idx_met}, met, 0.1*met, 'o', 'Color', colors(1, :))
        hold on
        % Model
        t_start   = data_s.time_mets{idx_met}(1)*30 + data_s.T1;
        t_end     = data_s.time_mets{idx_met}(end)*30 + data_s.T1;
        temps_loc = t_start:0.1:t_end;
        if isKey(param_dict, '$\tau$') && (param_dict('$\tau$') > 0)
            %------------------------------------------------------------------
            % Dormancy
            %------------------------------------------------------------------
            tau       = param_dict('$\tau$');
            V         = model_s.growth_model_met(...
                param_all(model_s.param_s.mets_growth_idx),...
                temps_loc - (birth_times(idx_met) + tau),... % birth_times are ordered in the same way as mets kinetics
                model_s.S0);
            if (birth_times(idx_met) + tau) > temps_loc(1)
                plot([(temps_loc(1) - shift_x)/scale, (birth_times(idx_met) + tau - shift_x)/scale],...
                    [model_s.unit_met_model2data(1), model_s.unit_met_model2data(1)], 'r')
                temps_loc = (birth_times(idx_met) + tau):0.1:temps_loc(end);
            end
            hold on
        else
            %------------------------------------------------------------------
            % No dormancy
            %------------------------------------------------------------------
            V         = model_s.growth_model_met(...
                param_all(model_s.param_s.mets_growth_idx),...
                temps_loc - birth_times(idx_met),... % birth_times are ordered in the same way as mets kinetics
                model_s.S0);
        end
        figure(3)
        plot((temps_loc - shift_x)/scale, model_s.unit_met_model2data(V),...
            'color', colors(1, :))
        figure(4)
        plot((temps_loc - shift_x)/scale, model_s.unit_met_model2data(V),...
            'color', colors(1, :))
        hold off
        set(gca, 'XLim', [0, t_max])
        xlabel('Months post-diagnosis')
        ylabel('Metastases diameter')
        set_fonts_lines(gca)
        export_fig([folder_kinetics '/met_' num2str(idx_met) '.pdf']);
    end
    figure(3)
    hold off
    set(gca, 'XLim', [0, t_max])
    xlabel('Months post-diagnosis')
    ylabel('Metastases diameter')
    set_fonts_lines(gca)
    export_fig([folder_kinetics '/all_mets.pdf']);
end
function D = cell2diam(N)
% Converts number of cells into diameter in mm
D = vol2diam(N*1e-6);
function D = vol2diam(V)
% Converts volume into diameter assuming spherical shape
D = (6*V/pi).^(1/3);
function N = diam2cell(D)
% Converts diameter in mm into number of cells
N = 1/6*pi*D.^3*1e6;
