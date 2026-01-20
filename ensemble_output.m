clear; clc;

%% Settings
num_runs = 4; % num runs / seeds 
timetotal = 25;    % how many output files were written
base_dir = 'C:\VRES_2025\PhysiCell-1.14.2\output_ensemble';

%automating the loading of the outputs
A1 = 'output0000000';
A2 = 'output000000';
A3 = 'output00000'; 
A4 = 'output0000'; 
B = '.xml';


T = timetotal - 1;

% storage
Tcell_all = nan(T, num_runs);
Rod_all   = nan(T, num_runs);
time_all  = nan(T, num_runs);
ActiveT_all   = nan(T, num_runs);
InactiveT_all = nan(T, num_runs);

for r = 1:num_runs
    runDir = fullfile(base_dir, sprintf('run_%03d', r), 'output');

    for tcount = 2:timetotal
        if tcount < 11
            K = [A1 num2str(tcount-1,'%d') B];
        elseif tcount < 101
            K = [A2 num2str(tcount-1,'%d') B];
        elseif tcount < 1001
            K = [A3 num2str(tcount-1,'%d') B];
        else
            K = [A4 num2str(tcount-1,'%d') B];
        end
        MCDS = read_MultiCellDS_xml(K, runDir);

        % count cells by type
        types = MCDS.discrete_cells.metadata.type;

        T_c = (types == 0);
        R_c = (types == 1);

        Tcell_all(tcount-1, r) = sum(T_c);  
        Rod_all(tcount-1, r)   = sum(R_c);  

        time_all(tcount-1, r) = MCDS.metadata.current_time;


        % activation status for t cells
        state = MCDS.discrete_cells.custom.state;
        t_states = state(T_c);     

        ActiveT_all(tcount-1, r)   = sum(t_states > 0.5);
        InactiveT_all(tcount-1, r) = sum(t_states <= 0.5);
       
    end
end

time = time_all(:,1);

%% calculate mean and standard deviation
T_mean = mean(Tcell_all, 2, 'omitnan');
T_std  = std(Tcell_all, 0, 2, 'omitnan');
R_mean = mean(Rod_all, 2, 'omitnan');
R_std  = std(Rod_all, 0, 2, 'omitnan');
A_mean = mean(ActiveT_all, 2, 'omitnan');
A_std  = std(ActiveT_all, 0, 2, 'omitnan');
I_mean = mean(InactiveT_all, 2, 'omitnan');
I_std  = std(InactiveT_all, 0, 2, 'omitnan');

figure; hold on; grid on;


%%fill(x, y, color) draws a filled polygon
fill([time; flipud(time)], ...
     [T_mean-T_std; flipud(T_mean+T_std)], ...
     [0.3 0.6 1.0], 'FaceAlpha', 0.25, 'EdgeColor', 'none');

plot(time, T_mean, 'LineWidth', 2);

xlabel('Time (min)');
ylabel('T cell count');
title('T cell population: mean \pm std across runs');
set(gca,'FontSize',14);

%% Plot Rod cells (mean ± std)
figure; hold on; grid on;

fill([time; flipud(time)], ...
     [R_mean-R_std; flipud(R_mean+R_std)], ...
     [0.3 0.8 0.3], 'FaceAlpha', 0.25, 'EdgeColor', 'none');

plot(time, R_mean, 'LineWidth', 2);

xlabel('Time (min)');
ylabel('Rod cell count');
title('Rod population: mean \pm std across runs');
set(gca,'FontSize',14);



figure; hold on; grid on;

hAfill = fill([time; flipud(time)], ...
              [A_mean-A_std; flipud(A_mean+A_std)], ...
              [1.0 0.5 0.2], 'FaceAlpha', 0.20, 'EdgeColor', 'none');

hAline = plot(time, A_mean, 'Color', [1.0 0.5 0.2], 'LineWidth', 2);

hIfill = fill([time; flipud(time)], ...
              [I_mean-I_std; flipud(I_mean+I_std)], ...
              [0.2 0.5 1.0], 'FaceAlpha', 0.20, 'EdgeColor', 'none');

hIline = plot(time, I_mean, 'Color', [0.2 0.5 1.0], 'LineWidth', 2);

xlabel('Time (min)');
ylabel('T cell count');
title('Active vs Inactive T cells: mean \pm std');

legend([hAfill hAline hIfill hIline], ...
       {'Active \pm std','Active mean','Inactive \pm std','Inactive mean'}, ...
       'Location','best');

set(gca,'FontSize',14);

