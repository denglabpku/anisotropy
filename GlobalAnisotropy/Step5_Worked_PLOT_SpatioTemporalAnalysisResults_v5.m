%   PLOT_SpatioTemporalAnalysisResults.m
%   Anders Sejr Hansen, October 2017
%   Modified by Zuhui 2020/06/16 for publication

clearvars -except ParentPath SampleNamePrefix resultDir; clc; close all;

%   DESCRIPTION
%   This script will load a structured array with all the data from
%   SpatioTemporal angle analysis fpr several cell lines/samples and then
%   do a number of hopefully informative plots

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
%   1. Delete some figure, only keep angular polar histogram, angle histogram, mean
%   displacement vs fold change, time vs fold change, 1st/2nd displacement
%   heatmap.
%   2. Correct an error in weighted mean in mean displacement vs fold change
%   plot. Should use ##NumAnglesAtCloestMean## to weight Angle# at each
%   DistIter each dt.

%%
DataSet = 1; % number of final mat of SpatioTemporal analysis of anisotropy
c_map = flip(hot,1); % Define colormap of imagesc plot
c_map = c_map(20:200,:);
resultDir = '/home/dell/Documents/METHOD/denglabpku-repository/anisotropy/ExampleData/AngularAnalysis/GlobalAnisotropy/';

% Load the appropriate dataset
if DataSet == 1
    load([resultDir, filesep, 'SpatioTemporalAnalysis', filesep, 'PA646_SpatioTemporalAnalysis_2022-05-18_17_03_45.mat']);
    plot_path = [resultDir, filesep, 'QC_Plots', filesep];
    if ~exist(plot_path,'dir')
        mkdir(plot_path)
    end
    SaveAppend = 'HMMfirst';
end

% Rename loaded sample name
FinalResults(1,1).SampleName = 'U2OS FOXA2-Halo';

%%
% Loop over all the samples:
for SampleIter = 1:length(FinalResults)
    tic;
    close; close;
    
    disp(['===plotting figure number 1 of ', num2str(SampleIter), '-th sample in total ', num2str(length(FinalResults)), ' result===']);
    figure1 = figure('Visible', 'on', 'position',[10 100 850 1000]); %[x y width height] 
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% SUB-PLOT 1 to 5: 
    %%%%%%  Angle Polar Distribution for 133 Hz
    %%%%%%  Angle histogram for 133 Hz
    %%%%%%  time vs fold change
    %%%%%%  mean displacement vs fold change
    %%%%%%  1st/2nd displacement fold change heatmap at 133 Hz
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% SUB-PLOT 1: Angle Polar Histogram for 133 Hz
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% 133 Hz plot: TimeIter = 1
    TimeIter = 1;
    subplot(3,2,1);
    polarplot(FinalResults(1,SampleIter).Alla{TimeIter},FinalResults(1,SampleIter).AllVals{TimeIter},...
        'Color', 'black');
    ax = gca;
    ax.RLim = [0 0.04];
    ax.RGrid = false;
    ax.RAxis.Visible = false;
    title([FinalResults(1,SampleIter).SampleName],'FontName', 'Helvetica', 'Interpreter', 'none');
    ax.FontSize = 14;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% SUB-PLOT 2: Angle histogram for 133 Hz
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
    
    %%% 133 Hz plot: TimeIter = 1
    TimeIter = 1;
    subplot(3, 2, 2);
    mean_FWHM_std = std(FinalResults(1,SampleIter).jack_FWHM{1, TimeIter}); %std of FWHM at 133 Hz
    hold on;
    plot(x_theta(x_index_vector), FinalResults(1,SampleIter).normPDF{1,TimeIter}(y_index_vector), '-', 'LineWidth', 2, 'Color', 'r');
    title([FinalResults(1,SampleIter).SampleName],...
                'FontName', 'Helvetica', 'Interpreter', 'none');
    text(0.6, 0.04, ['FWHM=' num2str(round(FinalResults(1,SampleIter).mean_FWHM(1,TimeIter),1)) '\pm' num2str(round(mean_FWHM_std ,1)) '\circ'],'FontSize',14);
    xlabel('\theta / 2\pi', 'FontName', 'Helvetica');
    ylabel('Probability', 'FontName', 'Helvetica'); 
    set(gca,'FontSize',14);
    hold off;
    
    
    
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
    subplot(3,2,3);
    hold on;
    plot(time_vector, FinalResults(1,SampleIter).mean_f_180_0, 'k--', 'LineWidth', 1);
    errorbar(time_vector, FinalResults(1,SampleIter).mean_f_180_0, std_f_180_0, 'ko', 'MarkerSize', 6, 'MarkerFaceColor',[237/255, 28/255, 36/255]);
    plot([0 1.02*max(time_vector)], [1 1], 'k--', 'LineWidth', 1);
    title([FinalResults(1,SampleIter).SampleName], 'FontName', 'Helvetica', 'Interpreter', 'none');
    ylabel('$$\mathbf{Fold(\frac{180\pm30^{\circ}}{0\pm30^{\circ}})}$$','Interpreter','latex');
    xlabel('lag time (ms)',  'FontName', 'Helvetica');
%     max_y = 4;
%     if mean(FinalResults(1,SampleIter).mean_f_180_0) > max_y
%         max_y = 6.1;
%     end
%     axis([0 1.02*max(time_vector) 0 max_y]);
    set(gca,'FontSize',14);
    hold off;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% SUB-PLOT 4: mean displacement vs fold change %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Plot a given anisotropy metric over space for all frame rates
    % diff_to_add to plot the y value in the middle of each bin
    diff_to_add = .5*(input_struct.MovingThreshold(2)-input_struct.MovingThreshold(1));
    
    %%% PLOT THE f(180/0) VS. SPACE
    subplot(3,2,4);
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
    title([FinalResults(1,SampleIter).SampleName],'FontName', 'Helvetica', 'Interpreter', 'none');
    ylabel('$$\mathbf{Fold(\frac{180\pm30^{\circ}}{0\pm30^{\circ}})}$$','Interpreter','latex');
    xlabel('mean displacement (nm)', 'FontName', 'Helvetica');
    max_y = 4.5;
    if mean(FinalResults(1,SampleIter).mean_f_180_0) > max_y
        max_y = 6.1;
    end
    axis([0 1000*(diff_to_add+1.03*max(input_struct.MovingThreshold(1,end))) 0.9 max_y]);
    set(gca,'FontSize',14);
    hold off;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SUB-PLOT 5: 1st/2nd displacement fold change heatmap at 133 Hz
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Plot a given anisotropy metric over space for all frame rates
    diff_to_add = .5*(input_struct.MovingThreshold(2)-input_struct.MovingThreshold(1));
    
    % Plot a heatmap of the fold(180/0)
    subplot(3,2,5);
    imagesc(1000*(input_struct.MovingThreshold), 1000*(input_struct.MovingThreshold), FinalResults(1,SampleIter).f_180_0_JumpMatrix);
    xlim([1000*(min(input_struct.MovingThreshold)-diff_to_add) 1000*(max(input_struct.MovingThreshold)+diff_to_add)]);
    ylim([1000*(min(input_struct.MovingThreshold)-diff_to_add) 1000*(max(input_struct.MovingThreshold)+diff_to_add)]);
    %set(gca,'ytick',frame_rate_vector,'yticklabel',legend_info);
    xlabel('Second displacement (nm)',  'FontName', 'Helvetica');
    ylabel('First displacement (nm)', 'FontName', 'Helvetica');
    title([FinalResults(1,SampleIter).SampleName], 'FontName', 'Helvetica', 'Interpreter', 'none');
    caxis([0.9 3.1]);
    colormap(c_map);
    c = colorbar('Location','eastoutside');
    c.Label.String = 'f(180/0)';
    c.Label.Rotation = 0;
    c.Label.Position = [0.6935 3.3673];
    set(gca,'FontSize',14);
    
    % Done with all the plotting: now save a PDF:
    set(figure1,'Units','Inches');
    pos = get(figure1,'Position');
    set(figure1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(figure1, [plot_path, FinalResults(1,SampleIter).SampleName, '_', SaveAppend, '_Plot1.pdf'], '-dpdf','-r0');
    toc;
end
close all;
