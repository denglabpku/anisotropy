
% This is a script to draw a randomly selected bound-unbound mixed
% trajectories from a given sample. The bound segments of a trajectory
% are color coded so that blue to red represents the start to
% end. The black segments represent bound segments. All charts have the
% same width as defined.

% Zuhui Wang 2021/09/17


%% HMM Processed CellTrack Input
clear; clc; close all;

Aniso_path = 'C:\Users\Zuhui\OneDrive - 北京大学生物医学前沿创新中心\Documents\MATLAB\SPT_Analysis\anisotropy-denglabpku\ExampleData\AngularAnalysis\LocalAnisotropy\';
Aniso_mat = 'anisotropy_maxFreeAfBfbound_FOXA2related_2022-05-13_21_05_09.mat';
load([Aniso_path Aniso_mat],'FinalResults','inputStruct');

fig_savepath = 'C:\Users\Zuhui\OneDrive - 北京大学生物医学前沿创新中心\Documents\MATLAB\SPT_Analysis\anisotropy-denglabpku\ExampleData\AngularAnalysis\LocalAnisotropy\Figures\';
DrawParams.ContinueFreeSegments = 4; % how many free segments before/after bound segments to plot
DrawParams.PlotWidth = 2; % unit: um
DrawParams.PlotNumber = 16; % tptal number of visualized randomly selected trajectories
DrawParams.MinMinJumpThres = 0.1;
DrawParams.MaxJump = 0.55 + 0.3; 
DrawParams.CircleRadius = 0.2; % circle radius around the end locs of bound segments, units: um.

%%      
for sampleIter =  [1 6]

    % load HMM classification data
    load([inputStruct(sampleIter).dataPath inputStruct(sampleIter).dataName]);
    
    % Determine which track are mixed with bound and unbound segments
    Temp_BoUnboMixTrackIdx = cellfun(@CheckUnique,CellTrackViterbiClass);
    %%% find track that are bound-free mixed
    BoUnboMixTrackIdx = find(Temp_BoUnboMixTrackIdx);
    BoUnbo_CellTrackViterbiClass = CellTrackViterbiClass(BoUnboMixTrackIdx);
    BoUnbo_CellTracks = CellTracks(BoUnboMixTrackIdx);
    fig_name = inputStruct(sampleIter).dataName(1:strfind(inputStruct(sampleIter).dataName,'Hz')+1);
    
    %% ONE BOUND FOLLOWED BY FREE SEGMENTS
    
    [FreAfBo_Single_f_180_0,PlotAf_CellTracks_idx,PlotAf_segment_idx] = SelectDrawIdx_FreeSegAfBo(sampleIter,FinalResults,DrawParams,BoUnbo_CellTracks);

    % Plot one bound segments followed by 5 free segments segment, colorcoded trajectory in single plot
    % Black line indicate the segment that are free
    % Blue to red indicate the appearance sequence of segments
    close all;
    DrawParams.PlotMode = 1; % Draw Bound followed by Free segments
    fig = figure('Visible','off');
    t = tiledlayout(4,4,'TileSpacing','compact','Padding','compact');
    for i = 1:DrawParams.PlotNumber     
        nexttile
        PlotColorCodeTrackSegment_v2(PlotAf_CellTracks_idx(i),PlotAf_segment_idx{i},BoUnbo_CellTracks,DrawParams)
    end

    fig.Units = 'inches';
    fig.Position = [4.3125 1.9479 9.4688 7.8750];
    t.Title.String = {inputStruct(sampleIter).dataName(1:end-4),['f_180_0:' num2str(FreAfBo_Single_f_180_0)]} ;
    t.Title.Interpreter = 'none';

    print([fig_savepath fig_name '_' num2str(DrawParams.ContinueFreeSegments) 'FreeSegAfterBo_' char(datetime('now','format','dd-MM-yyy_HH_mm_ss'))],'-dpdf','-r0','-bestfit')
    
    clear Plot*
    %% FREE SEGMENTS FOLLOWED BY ONE BOUND 
    [FreBfBo_Single_f_180_0, PlotBf_CellTracks_idx,PlotBf_segment_idx] = SelectDrawIdx_FreeSegBfBo(sampleIter,FinalResults,DrawParams,BoUnbo_CellTracks);
    
    % Plot 5 free segments segments followed by one bound, colorcoded trajectory in single plot
    % Black line indicate the segment that are free
    % Blue to red indicate the appearance sequence of segments
    close all;
    DrawParams.PlotMode = 2; % Draw free segments followed by Bound  
    fig = figure('Visible','off');
    t = tiledlayout(4,4,'TileSpacing','compact','Padding','compact');
    for i = 1:DrawParams.PlotNumber     
        nexttile
        PlotColorCodeTrackSegment_v2(PlotBf_CellTracks_idx(i),PlotBf_segment_idx{i},BoUnbo_CellTracks,DrawParams)
    end

    fig.Units = 'inches';
    fig.Position = [4.3125 1.9479 9.4688 7.8750];
    t.Title.String = {inputStruct(sampleIter).dataName(1:end-4),['f_180_0:' num2str(FreBfBo_Single_f_180_0)]};
    t.Title.Interpreter = 'none';

    print([fig_savepath fig_name '_' num2str(DrawParams.ContinueFreeSegments) 'FreeSegBeforeBo_' char(datetime('now','format','yyyy-MM-dd_HH_mm_ss'))],'-dpdf','-r0','-bestfit')
    clear Plot*
    
    clear CellTracks CellTrackViterbiClass
end

%% AUXILIARY FUNCTION %%

function ifBoUnboMix = CheckUnique(x)
    State = unique(x);
     
    if sum(State) == 3 % Bound/unbound segment mixed 
        ifBoUnboMix = true;
    else
        ifBoUnboMix = false;
    end
end

function FreeSeg_count = CountFreeSeg(x)
    FreeSeg_count = sum(x == 2);
end

function BoundSeg_count = CountBoundSeg(x)
    BoundSeg_count = sum(x == 1);
end





