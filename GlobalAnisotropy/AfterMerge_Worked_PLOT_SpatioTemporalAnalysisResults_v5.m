%   PLOT_SpatioTemporalAnalysisResults.m
%   Anders Sejr Hansen, October 2017
%   Modified by Zuhui 2021/05/10

clear; clc; close all;

%   DESCRIPTION This script will load a structured array with all the data
%   merged FinalResults from multiple cell lines and then do a number of
%   hopefully informative plots

%   UPDATE - v2
%   When plot anisotropy metrics over space, average over all timepoints.

%   UPDATE - v4
%   When analyzing displacements over lengths, you will have much fewer
%   angles for the long lag times. So weigh the angle metrics according to
%   the number of angles for each time point
%   Also plot metrics for MIN instead of MEAN displacement
%   Finally plot the full spatiotemporal matrix, where the first dimension
%   is the first jump and the second dimension is the second jump

%   UPDATE - v5
%   Delete some figure, only keep angular polar histogram, angle histogram, mean
%   displacement vs fold change, time vs fold change, 1st/2nd displacement
%   heatmap.
%   Correct an error in weighted mean in mean displacement vs fold change
%   plot. Should use ##NumAnglesAtCloestMean## to weight Angle# at each
%   DistIter each dt.

%% Merge FinalResults from individual cell line...
% After get the FinalResults from MultipleFrameRate script, manually
% horizontally concatenate FinalResults

clc, clear, close all;
DataSet = 1; % number of final mat of SpatioTemporal analysis of anisotropy
c_map = flip(hot,1); % Define colormap of imagesc plot
c_map = c_map(20:200,:);
resultDir = '/mnt/f58069a5-1cf3-43b8-bb9b-ea74327327c9/WZH-DataCenter/PROCESS-SPT/AnisotropyPlot_XLONEFOXA2related';


% Rename loaded sample name
% Load the appropriate dataset
if DataSet == 1
    load(fullfile(resultDir, 'Merged2TimePoints_PA646_SpatioTemporalAnalysis_2022-04-14_15_41_37.mat'));
    plot_path = [resultDir, filesep, '20220414_Merged2TimePoints_2022-04-14_15_41_37', filesep];
    if ~exist(plot_path,'dir')
        mkdir(plot_path);
    end
    SaveAppend = '2TimePoints_2022-04-14_15_41_37';
end

FinalResults = Merge_FinalResults;
FinalResults(1,1).SampleName = 'U2OS FOXA2-Halo';
FinalResults(1,2).SampleName = 'U2OS \DeltaNTD-Halo';
FinalResults(1,3).SampleName = 'U2OS \DeltaDBD-Halo';
FinalResults(1,4).SampleName = 'U2OS \DeltaCTD-Halo';
FinalResults(1,5).SampleName = 'U2OS DBD-Halo';

input_struct.MovingThreshold = Merge_input_struct(1).MovingThreshold;

%% Loop over all the samples:
for SampleIter = 1:length(FinalResults)
    tic;
    close; close;
    
    disp(['===plotting figure number 1 of ', num2str(SampleIter), '-th sample in total ', num2str(length(FinalResults)), ' result===']);
    %figure1 = figure('position',[10 100 600 1000]); %[x y width height] 
    figure1 = figure('Visible', 'off'); %[x y width height] 
    tiledlayout(3,3,"TileSpacing","tight");
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% SUB-PLOT 1 to 5: 
    %%%%%%  Angle Polar Distribution for 222 Hz
    %%%%%%  Angle histogram for 222 Hz
    %%%%%%  time vs fold change
    %%%%%%  mean displacement vs fold change
    %%%%%%  1st/2nd displacement fold change heatmap at 222 Hz
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% SUB-PLOT 1 to 3: ANGLE DISTRIBUTION FOR 223, 133 and 74 Hz
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 223 Hz plot: TimeIter = 1
    % 133 Hz plot: TimeIter = 2
    % 74 Hz plot: TimeIter = 3
    for TimeIter = 1:3
        nexttile;
        %Option 1:
        %Max limit:
    %     MaxLim = max([0.04 1.01*max(FinalResults(1,SampleIter).AllVals{TimeIter})]);
    %     h = polar(0,MaxLim);
    %     delete(h);
    % %     r = groot;
    % %     r.ShowHiddenHandles = 'on';
    % %     findobj(gca, 'Type','text')... to invisiblize RGrid
    %     set(gca, 'Nextplot','add')
    %     %# draw patches instead of lines: polar(t,r)
    %     [x,y] = pol2cart(FinalResults(1,SampleIter).Alla{TimeIter},FinalResults(1,SampleIter).AllVals{TimeIter});
    %     h = patch(reshape(x,4,[]), reshape(y,4,[]), [237/255, 28/255, 36/255]);
    %     title([FinalResults(1,SampleIter).SampleName],...
    %             'FontSize', 8, 'FontName', 'Helvetica', 'Interpreter', 'none');
    %     xlabel('angle (degrees)', 'FontSize',9, 'FontName', 'Helvetica');

        %Option 2:
        polarplot(FinalResults(1,SampleIter).Alla{TimeIter},FinalResults(1,SampleIter).AllVals{TimeIter},...
            'Color', 'black');
        ax = gca;
        ax.RLim = [0 0.04];
        ax.RGrid = false;
        ax.RAxis.Visible = false;
        title({[FinalResults(1,SampleIter).SampleName, FinalResults(1,SampleIter).FrameRates(1,TimeIter).FrameRate, ': angles'];...
                ['AC = ', num2str(FinalResults(1,SampleIter).mean_AC(1,TimeIter)), ...
                    ' +/- ', num2str(std(FinalResults(1,SampleIter).jack_AC{TimeIter}))];...
                ['f(180/0) = ', num2str(FinalResults(1,SampleIter).mean_f_180_0(1,TimeIter)), ...
                    ' +/- ', num2str(std(FinalResults(1,SampleIter).jack_f_180_0{TimeIter}))];}...
                ,'FontSize', 8, 'FontName', 'Helvetica');
        ax.FontSize = 8;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% SUB-PLOT 2: Angle histogram for 223, 133 and 74 Hz %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % first generate the required index_vectors:
    normPDF = FinalResults(1,SampleIter).normPDF{1};
    % x-vector for displaying the angles
    x_index_vector = [1;];
    for i=2:length(normPDF)
        x_index_vector = [x_index_vector, i, i];
    end
    x_index_vector = [x_index_vector, length(normPDF)+1];
    y_index_vector = [];
    for i=1:length(normPDF)
        y_index_vector = [y_index_vector, i, i];
    end
    x_theta = 0:1/length(normPDF):1;
    
    for TimeIter = 1:3
        nexttile;
        mean_FWHM_std = std(FinalResults(1,SampleIter).jack_FWHM{1, TimeIter}); %std of FWHM at 222 Hz
        hold on;
        plot(x_theta(x_index_vector), FinalResults(1,SampleIter).normPDF{1,TimeIter}(y_index_vector), '-', 'LineWidth', 2, 'Color', 'r');
        title({['AC = ', num2str(FinalResults(1,SampleIter).mean_AC(1,TimeIter)), ...
                    ' +/- ', num2str(std(FinalResults(1,SampleIter).jack_AC{TimeIter}))];...
                ['f(180/0) = ', num2str(FinalResults(1,SampleIter).mean_f_180_0(1,TimeIter)), ...
                    ' +/- ', num2str(std(FinalResults(1,SampleIter).jack_f_180_0{TimeIter}))];}...
                ,'FontSize', 8, 'FontName', 'Helvetica');
        xlabel('\theta / 2\pi', 'FontName', 'Helvetica');
        ylabel('Probability', 'FontName', 'Helvetica'); 
        axis([0 1 0.02 0.05]);
        set(gca,'FontSize',8);
        hold off;
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% SUB-PLOT 3: time vs fold change
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Plot f(180/0) as a function of time with error bars
    % First need to get the standard deviations
%     std_Amp = zeros(1,length(FinalResults(1,SampleIter).mean_Amp));
    std_f_180_0 = zeros(1,length(FinalResults(1,SampleIter).mean_f_180_0));
%     std_FWHM = zeros(1,length(FinalResults(1,SampleIter).mean_FWHM));
    time_vector = zeros(1,length(FinalResults(1,SampleIter).FrameRates));
    for TimeIter = 1:length(std_f_180_0)
%         std_Amp(1,TimeIter) = std(FinalResults(1,SampleIter).jack_Amp{TimeIter});
        std_f_180_0(1,TimeIter) = std(FinalResults(1,SampleIter).jack_f_180_0{TimeIter});
%         std_FWHM(1,TimeIter) = std(FinalResults(1,SampleIter).jack_FWHM{TimeIter});
        time_vector(1,TimeIter) = 1000*FinalResults(1,SampleIter).FrameRates(1,TimeIter).LagTime;
    end
     
    % Plot f(180/0) vs. time
    nexttile;
    hold on;
    plot(time_vector, FinalResults(1,SampleIter).mean_f_180_0, 'k--', 'LineWidth', 1);
    errorbar(time_vector, FinalResults(1,SampleIter).mean_f_180_0, std_f_180_0, 'ko', 'MarkerSize', 6, 'MarkerFaceColor',[237/255, 28/255, 36/255]);
    plot([0 1.02*max(time_vector)], [1 1], 'k--', 'LineWidth', 1);
    title([FinalResults(1,SampleIter).SampleName], 'FontName', 'Helvetica');
    ylabel('$$\mathbf{Fold(\frac{180\pm30^{\circ}}{0\pm30^{\circ}})}$$','Interpreter','latex');
    xlabel('lag time (ms)',  'FontName', 'Helvetica');
    max_y = 5;
    if mean(FinalResults(1,SampleIter).mean_f_180_0) > max_y
        max_y = 6.1;
    end
    axis([0 1.02*max(time_vector) 0 max_y]);
    set(gca,'FontSize',14);
    hold off;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% SUB-PLOT 4: mean displacement vs fold change %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Plot a given anisotropy metric over space for all frame rates
    % diff_to_add to plot the y value in the middle of each bin
    diff_to_add = .5*(input_struct.MovingThreshold(2)-input_struct.MovingThreshold(1));
    
    %%% PLOT THE f(180/0) VS. SPACE
    nexttile;
    hold on;
    %%%%%%% CALCULATE USING WEIGHTS BASED ON NUMBER OF ANGLES %%%%%%%
    
    % Description of calculation:
    % We want to calculate the f(180/0) as a function of the mean
    % displacement, averaged over all the different deltaT. 
    % From the processed calculations we have a matrix: 
    % FinalResults.f_180_0_AtClosestMean. 
    % This is an TxL matrix, where:
    % There are T columns for the T timepoints: 223 Hz, 134 Hz, ... 9.2 Hz
    % There are L rows for the L displacement bins (contained in the "input_struct.MovingThreshold" vector
    % So here we want to take the mean anisotropy for a given mean
    % displacement length bin L and average over all the timepoints T. 
    % But here there are 2 issues to consider:
    % - empty elements (e.g. if there was not enough angles to accurately
    % calculate the f(180/0) it will be left empty or if it was restricted
    % (e.g. for 223 Hz data the f(180/0) is not calculated for
    % displacements > 450 nm and therefore empty at these distances)
    % - the number of angles are different. For example, if we calculate
    % the f(180/0) for 200 nm bins and the 223 Hz dataset contributes 50k
    % angles, but the 9.2 Hz dataset contributes 5k angles, we should weigh
    % them accordingly when we take the mean (in this examples, the 223 Hz
    % dataset should count 10x of the 9.2 Hz dataset). 
    % The code snippet below therefore does this:
    % It first removes any empty (NaN) values. It then calculates a
    % weighing vector "curr_NumAngles", which sums to 1. Finally, it
    % re-scales the f(180/0) value according to this. 
    
    dist_mean = zeros(1,length(input_struct.MovingThreshold));
    dist_std = zeros(1,length(input_struct.MovingThreshold));
    for DistIter = 1:length(dist_mean)
        curr_vals = FinalResults(1,SampleIter).f_180_0_AtClosestMean(:,DistIter);
        % should use ##NumAnglesAtCloestMean## to weight Angle# at each DistIter each dt
        curr_NumAngles = FinalResults(1,SampleIter).NumAnglesAtClosestMean(:,DistIter);
%         curr_NumAngles = FinalResults(1,SampleIter).TotalMinMinAngles'; % convert to column vector
        % get rid of NaN's to be able to calculate mean/std
        ToRemove = isnan(curr_vals);
        curr_vals(ToRemove) = [];
        % remove the same elements in the TotalMinMinAngles counter:
        curr_NumAngles(ToRemove) = [];
        % normalize the curr_NumAngles:
        dist_mean(1,DistIter) = sum(curr_vals.*curr_NumAngles) / sum(curr_NumAngles);
        % SYNTAX of STD in MATLAB 2014b: std( value_vector, weight_vector );
        dist_std(1,DistIter) = std(curr_vals, curr_NumAngles)/ sqrt(length(curr_vals));
    end
    
    % now do the actual plotting:
    plot(1000*(input_struct.MovingThreshold+diff_to_add),dist_mean, 'k--','LineWidth',1);
    errorbar(1000*(input_struct.MovingThreshold+diff_to_add),dist_mean, dist_std, 'ko', 'MarkerSize', 6,  'MarkerFaceColor', [237/255, 28/255, 36/255]);
    plot([0 5000], [1 1], 'k--', 'LineWidth', 1);
    title([FinalResults(1,SampleIter).SampleName],'FontName', 'Helvetica');
    ylabel('$$\mathbf{Fold(\frac{180\pm30^{\circ}}{0\pm30^{\circ}})}$$','Interpreter','latex');
    xlabel('mean displacement (nm)', 'FontName', 'Helvetica');
    max_y = 5;
    if mean(FinalResults(1,SampleIter).mean_f_180_0) > max_y
        max_y = 6.1;
    end
    axis([0 1000*(diff_to_add+1.03*max(input_struct.MovingThreshold(1,end))) 0.9 max_y]);
    set(gca,'FontSize',14);
    hold off;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SUB-PLOT 5: 1st/2nd displacement fold change heatmap at 222 Hz
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Plot a given anisotropy metric over space for all frame rates
    diff_to_add = .5*(input_struct.MovingThreshold(2)-input_struct.MovingThreshold(1));
    
    % Plot a heatmap of the fold(180/0)
    nexttile;
    imagesc(1000*(input_struct.MovingThreshold), 1000*(input_struct.MovingThreshold), FinalResults(1,SampleIter).f_180_0_JumpMatrix);
    xlim([1000*(min(input_struct.MovingThreshold)-diff_to_add) 1000*(max(input_struct.MovingThreshold)+diff_to_add)]);
    ylim([1000*(min(input_struct.MovingThreshold)-diff_to_add) 1000*(max(input_struct.MovingThreshold)+diff_to_add)]);
    %set(gca,'ytick',frame_rate_vector,'yticklabel',legend_info);
    xlabel('Second displacement (nm)',  'FontName', 'Helvetica');
    ylabel('First displacement (nm)', 'FontName', 'Helvetica');
    title([FinalResults(1,SampleIter).SampleName], 'FontName', 'Helvetica');
    caxis([0.9 3.1]);
    colormap(c_map);
    c = colorbar('Location','eastoutside');
    c.Label.String = 'f(180/0)';
    c.Label.Rotation = 0;
    c.Label.Position = [0.6935 3.3673];
    set(gca,'FontSize',14);
    
    % Done with all the plotting: now save a PDF:
    set(figure1,'Units','Inches');
    pos = [3.1875   -0.1562   12.3646   10.4167];
    set(figure1,'PaperSize',[pos(3), pos(4)],'Position',pos);
    print(figure1, [plot_path, Merge_FinalResults(1,SampleIter).SampleName, SaveAppend, '_PlotSingle.pdf'], '-dpdf','-r0');
    toc;
end
close all;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% PROCEED TO THE 2ND FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('===plotting overlay figure===');
figure2 = figure('Visible', 'off'); %[x y width height]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% SUB-PLOT 1 to 3: 
%%%%%%  Overlaid Angle histogram for 222 Hz
%%%%%%  Overlaid time vs fold change
%%%%%%  Overlaid mean displacement vs fold change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sample color choices of overlaid plot
ColorChoice = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE'};
tiledlayout(2,3,"TileSpacing","tight");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% SUB-PLOT 1: Overlaid Angle histogram for 223, 133 and 74 Hz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for TimeIter = 1:3
    nexttile;
    hold on;
    iter = 1;
    legendIter = 1:length(FinalResults);
    for SampleIter = 1:length(FinalResults)
        % first generate the required index_vectors:
        normPDF = FinalResults(1,SampleIter).normPDF{1};
            % x-vector for displaying the angles
            x_index_vector = [1;];
            for i=2:length(normPDF)
                x_index_vector = [x_index_vector, i, i];
            end
            x_index_vector = [x_index_vector, length(normPDF)+1];
            y_index_vector = [];
            for i=1:length(normPDF)
                y_index_vector = [y_index_vector, i, i];
            end
            x_theta = 0:1/length(normPDF):1;

        mean_FWHM_std = std(FinalResults(1,SampleIter).jack_FWHM{1, TimeIter}); %std of FWHM at 222 Hz
        plot(x_theta(x_index_vector), FinalResults(1,SampleIter).normPDF{1,TimeIter}(y_index_vector), '-', 'LineWidth', 2, 'Color', ColorChoice{iter});
        text(0.02, 0.047-iter*0.002, ['FWHM=' num2str(round(FinalResults(1,SampleIter).mean_FWHM(1,TimeIter),1)) '\pm' num2str(round(mean_FWHM_std ,1)) '\circ'],'FontSize',14,'Color', ColorChoice{iter});
        iter = iter +1;
    end
    clear iter;
    axis([0 1 0.02 0.05]);
    title([FinalResults(1,1).FrameRates(1,TimeIter).FrameRate, ': histogram']);
    xlabel('\theta / 2\pi', 'FontName', 'Helvetica');
    ylabel('Probability', 'FontName', 'Helvetica');
    legend(FinalResults(1,legendIter).SampleName,'Box','off','FontSize',10);
    set(gca,'FontSize',14);
    hold off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% SUB-PLOT 2: Overlaid time vs fold change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile;
hold on;
iter = 1;
for SampleIter = 1:length(FinalResults) 
% Plot f(180/0) as a function of time with error bars
    % First need to get the standard deviations
%     std_Amp = zeros(1,length(FinalResults(1,SampleIter).mean_Amp));
    std_f_180_0 = zeros(1,length(FinalResults(1,SampleIter).mean_f_180_0));
%     std_FWHM = zeros(1,length(FinalResults(1,SampleIter).mean_FWHM));
    time_vector = zeros(1,length(FinalResults(1,SampleIter).FrameRates));
    for TimeIter = 1:length(std_f_180_0)
%         std_Amp(1,TimeIter) = std(FinalResults(1,SampleIter).jack_Amp{TimeIter});
        std_f_180_0(1,TimeIter) = std(FinalResults(1,SampleIter).jack_f_180_0{TimeIter});
%         std_FWHM(1,TimeIter) = std(FinalResults(1,SampleIter).jack_FWHM{TimeIter});
        time_vector(1,TimeIter) = 1000*FinalResults(1,SampleIter).FrameRates(1,TimeIter).LagTime;
    end
     
    % Plot f(180/0) vs. time
    plot(time_vector, FinalResults(1,SampleIter).mean_f_180_0, 'k--', 'LineWidth', 1,'HandleVisibility','off');
    errorbar(time_vector, FinalResults(1,SampleIter).mean_f_180_0, std_f_180_0, 'ko',...
        'LineWidth', 1, 'MarkerSize', 6, 'MarkerFaceColor',ColorChoice{iter});
    max_y = 4;
    if mean(FinalResults(1,SampleIter).mean_f_180_0) > max_y
        max_y = 5;
    end
    iter = iter +1;
end
clear iter;
axis([0 1.02*max(time_vector) 0 max_y]);
title('anisotropy: time dependency', 'FontName', 'Helvetica', 'Interpreter', 'none');
ylabel('$$\mathbf{Fold(\frac{180\pm30^{\circ}}{0\pm30^{\circ}})}$$','Interpreter','latex');
xlabel('lag time (ms)', 'FontName', 'Helvetica');
legend(FinalResults(1,legendIter).SampleName,'Box','off');
plot([0 1.02*max(time_vector)], [1 1], 'k--', 'LineWidth', 1,'HandleVisibility','off');
set(gca,'FontSize',14);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% SUB-PLOT 3: Overlaid mean displacement vs fold change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot a given anisotropy metric over space for all frame rates
diff_to_add = .5*(input_struct.MovingThreshold(2)-input_struct.MovingThreshold(1));
    
%%% PLOT THE f(180/0) VS. SPACE
nexttile;
%%%%%%% CALCULATE USING WEIGHTS BASED ON NUMBER OF ANGLES %%%%%%%

% Description of calculation:
% We want to calculate the f(180/0) as a function of the mean
% displacement, averaged over all the different deltaT. 
% From the processed calculations we have a matrix: 
% FinalResults.f_180_0_AtClosestMean. 
% This is an TxL matrix, where:
% There are T columns for the T timepoints: 223 Hz, 134 Hz, ... 9.2 Hz
% There are L rows for the L displacement bins (contained in the "input_struct.MovingThreshold" vector
% So here we want to take the mean anisotropy for a given mean
% displacement length bin L and average over all the timepoints T. 
% But here there are 2 issues to consider:
% - empty elements (e.g. if there was not enough angles to accurately
% calculate the f(180/0) it will be left empty or if it was restricted
% (e.g. for 223 Hz data the f(180/0) is not calculated for
% displacements > 450 nm and therefore empty at these distances)
% - the number of angles are different. For example, if we calculate
% the f(180/0) for 200 nm bins and the 223 Hz dataset contributes 50k
% angles, but the 9.2 Hz dataset contributes 5k angles, we should weigh
% them accordingly when we take the mean (in this examples, the 223 Hz
% dataset should count 10x of the 9.2 Hz dataset). 
% The code snippet below therefore does this:
% It first removes any empty (NaN) values. It then calculates a
% weighing vector "curr_NumAngles", which sums to 1. Finally, it
% re-scales the f(180/0) value according to this. 
hold on;
iter = 1;
for SampleIter = 1:length(FinalResults) 
    dist_mean = zeros(1,length(input_struct.MovingThreshold));
    dist_std = zeros(1,length(input_struct.MovingThreshold));
    for DistIter = 1:length(dist_mean)
        curr_vals = FinalResults(1,SampleIter).f_180_0_AtClosestMean(:,DistIter);
        % should use ##NumAnglesAtCloestMean## to weight Angle# at each DistIter each dt
        curr_NumAngles = FinalResults(1,SampleIter).NumAnglesAtClosestMean(:,DistIter);
        % curr_NumAngles = FinalResults(1,SampleIter).TotalMinMinAngles'; % convert to column vector
        % get rid of NaN's to be able to calculate mean/std
        ToRemove = isnan(curr_vals);
        curr_vals(ToRemove) = [];
        % remove the same elements in the TotalMinMinAngles counter:
        curr_NumAngles(ToRemove) = [];
        % normalize the curr_NumAngles:
        dist_mean(1,DistIter) = sum(curr_vals.*curr_NumAngles) / sum(curr_NumAngles);
        % SYNTAX of STD in MATLAB 2014b: std( value_vector, weight_vector );
        dist_std(1,DistIter) = std(curr_vals, curr_NumAngles)/ sqrt(length(curr_vals));
    end
    
    % now do the actual plotting:
    plot(1000*(input_struct.MovingThreshold(1:12)+diff_to_add),dist_mean(1:12), 'k--','LineWidth',1, 'HandleVisibility','off');
    errorbar(1000*(input_struct.MovingThreshold(1:12)+diff_to_add),dist_mean(1:12), dist_std(1:12), 'ko',...
        'MarkerSize', 6,  'MarkerFaceColor', ColorChoice{iter});
    iter = iter +1;
end
max_y = 4.5;
clear iter;
% axis([0 1000*(diff_to_add+1.03*max(input_struct.MovingThreshold(1,12))) 0 max_y]);
axis([0 700 0 max_y]);
title('anisotropy: space dependency', 'FontName', 'Helvetica', 'Interpreter', 'none');
ylabel('$$\mathbf{Fold(\frac{180\pm30^{\circ}}{0\pm30^{\circ}})}$$','Interpreter','latex');
xlabel('mean displacement (nm)', 'FontName', 'Helvetica');
legend(FinalResults(1,legendIter).SampleName,'Box','off');
plot([0 5000], [1 1], 'k--', 'LineWidth', 1, 'HandleVisibility','off');
set(gca,'FontSize',14);
hold off;

% Done with all the plotting: now save a PDF:
set(figure2,'Units','Inches');
pos = [0.7500    0.0104   20   10.3438];
set(figure2,'PaperSize',[pos(3), pos(4)],'Position',pos)
print(figure2, [plot_path, 'Overlay', '_', SaveAppend], '-dpdf','-r0');
disp('===Overlay plot export completed===');


