% This is a script to plot the local anisotropy at the given free segment
% lengths either before or after a bound segment.

% Zuhui Wang 2021/10/17

close all; clear; clc;
%% Plot local anisotropy
%% Plot the max_FreSeg_size vs f_180_0
input_path = 'C:\Users\Zuhui\OneDrive - 北京大学生物医学前沿创新中心\Documents\MATLAB\SPT_Analysis\anisotropy-denglabpku\ExampleData\AngularAnalysis\LocalAnisotropy\';
input_name = 'anisotropy_maxFreeAfBfbound_FOXA2related_2022-05-13_21_05_09.mat';
load([input_path input_name],'FinalResults','inputStruct','max_FreSeg_size');

fig_output = [input_path filesep 'Figures'];

colorChoice = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE'};

fig = figure('Visible','off','PaperOrientation','landscape');
tiledlayout(2,2,'TileSpacing','compact');

%==plot title=== anisotropy of free segments after one bound segment
nexttile
hold on
title('anisotropy of free segments after one bound segment')
xlabel('free segment lengths')
ylabel('$$\mathbf{Fold(\frac{180\pm30^{\circ}}{0\pm30^{\circ}})}$$','Interpreter','latex');
for SampleIter = 1:length(inputStruct)

    errorbar(max_FreSeg_size,[FinalResults(:,SampleIter).FreAfBo_mean_f_180_0],[FinalResults(:,SampleIter).FreAfBo_std_f_180_0],...
        'LineStyle','--','Color','k','LineWidth',1,...
        'Marker','o','MarkerSize',6,'MarkerFaceColor',colorChoice{SampleIter},'MarkerEdgeColor','k');
end
xlim([1 max(max_FreSeg_size)+1]);
ylim([1 4.5]);
xticks(1:max(max_FreSeg_size)+1);
hold off
set(gca,'FontSize',14);

%==plot title=== anisotropy of free segments before one bound segment
nexttile
hold on
title('anisotropy of free segments before one bound segment')
xlabel('free segment lengths')
ylabel('$$\mathbf{Fold(\frac{180\pm30^{\circ}}{0\pm30^{\circ}})}$$','Interpreter','latex');
for SampleIter = 1:length(inputStruct)
    errorbar(max_FreSeg_size,[FinalResults(:,SampleIter).FreBfBo_mean_f_180_0],[FinalResults(:,SampleIter).FreBfBo_std_f_180_0],...
        'LineStyle','--','Color','k','LineWidth',1,...
        'Marker','o','MarkerSize',6,'MarkerFaceColor',colorChoice{SampleIter},'MarkerEdgeColor','k');
end
xlim([1 max(max_FreSeg_size)+1]);
ylim([1 4.5]);
xticks(1:max(max_FreSeg_size)+1);
% text_legend = {inputStruct(:).sample_prefix};
% legend(text_legend,'Box','off');
set(gca,'FontSize',14);

nexttile
hold on
title('anisotropy of free segments after or before bound segments')
xlabel('free segment lengths')
ylabel('$$\mathbf{Fold(\frac{180\pm30^{\circ}}{0\pm30^{\circ}})}$$','Interpreter','latex');
for SampleIter = 1:length(inputStruct)
    errorbar(max_FreSeg_size,[FinalResults(:,SampleIter).FreABfBo_mean_f_180_0],[FinalResults(:,SampleIter).FreABfBo_std_f_180_0],...
        'LineStyle','--','Color','k','LineWidth',1,...
        'Marker','o','MarkerSize',6,'MarkerFaceColor',colorChoice{SampleIter},'MarkerEdgeColor','k');
end
xlim([1 max(max_FreSeg_size)+1]);
ylim([1 4.5]);
xticks(1:max(max_FreSeg_size)+1);
text_legend = {inputStruct(:).sample_prefix};
legend(text_legend,'Box','off');
set(gca,'FontSize',14);

fig.Units = 'inches';
fig.Position = [6,4,9,8];
fig.PaperType = 'A4';

fig_name = 'anisotropy_maxFreeAfBfbound_FOXA2related_globalAdd';
print(fullfile(fig_output, fig_name),'-dpdf','-r0'); % open svg in illustrator for further editing

%% Plot the polar histogram of angles from pooled global free segments
clear; close all;
input_path = 'C:\Users\Zuhui\OneDrive - 北京大学生物医学前沿创新中心\Documents\MATLAB\SPT_Analysis\anisotropy-denglabpku\ExampleData\AngularAnalysis\LocalAnisotropy\';
input_name = 'anisotropy_maxFreeAfBfbound_FOXA2related_2022-05-13_21_05_09.mat'; 
load([input_path input_name],'FinalResults','inputStruct','max_FreSeg_size');

fig_output = [input_path filesep 'Figures'];

for segIter = [1 2] % max_segment 2 and 3
    fig = figure('Visible','off');
    t = tiledlayout(3,2,'TileSpacing','compact');
    title(t,sprintf('Polar histogram (free segment length: %d)',FinalResults(segIter,1).max_FreSeg_size));
    for SampleIter = 1:6
        nexttile;
        polarplot(FinalResults(segIter,SampleIter).Alla,FinalResults(segIter,SampleIter).AllVals,...
                'Color', 'black');
        title([FinalResults(segIter,SampleIter).sample_prefix],...
            'FontName', 'Helvetica', 'Interpreter', 'none');
        ax = gca;
        ax.RLim = [0 0.04];
        ax.RGrid = false;
        ax.RAxis.Visible = false;
        ax.FontSize = 14;
    end
    fig_name = sprintf('polarhisto_FreeAfBfbound_FOXA2related_globalAdd_Plot%d',segIter);
    fig.Units = 'inches';
    fig.Position = [1.7129 0.6139 8.0000 9];
    fig.PaperType = 'A4';
    print(fullfile(fig_output, fig_name),'-dpdf','-r0'); % open svg in illustrator for further editing
end
fprintf('Writing figure complete \n');
