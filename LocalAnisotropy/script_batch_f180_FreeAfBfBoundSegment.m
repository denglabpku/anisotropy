% This is a script to find the free segment before and after a bound
% segment and calculate the corresponding f_180_0 in a multisample batch mode.

% Zuhui Wang 2021/09/16

%% Note:
% Anders' Anisotropy analysis code defined two threshold to reduce strang long trajectory artifact
% MinNumAngles = 5;
% MaxAsymAnglesFrac = 0.5; 

% We do not incorporate above two threshold in this script because we only
% focus on a fraction of free segment either after or before a bound
% segment. We also have a required free segment length
% (FreSeg_size), so the segment length are not that long on average.
% 'MaxAsymAnglesFrac' this fraction is not that reasonable if the segment
% length are very short, since if a molecule indead exhibit a strong
% anisotropic diffusion behavior, the first few track after released from a
% bound segment can exhibit highly asymetric angles.

% Explanation on max_FreSeg_size: this variable currently defined the exact length of
% Free segments used to calculate the anistropy etc. See the line 116 and
% line 160 of FreeSegmentAsymCalculation.m;


%% Part I Merge 133Hz and 74Hz HMM_first_QC_data classified data
clear;clc;close all;

dataPath_prefix = '/home/dell/Documents/METHOD/denglabpku-repository/anisotropy/ExampleData/';

inputStruct(1).dataPath = [dataPath_prefix '/AngularAnalysis/GlobalAnisotropy/'];
inputStruct(1).sample_prefix = 'U2OS_FOXA2-Halo_';

% inputStruct(2).dataPath = [dataPath_prefix '2021/20210413-PA646-U2OS-XLONE-FOXA2dNTD-Halo/AngularAnalysis/'];
% inputStruct(2).sample_prefix = 'U2OS_Xlone-FOXA2dNTD-Halo_';
% 
% inputStruct(3).dataPath = [dataPath_prefix '2021/20210506-PA646-U2OS-XLONE-FOXA2dDBD-Halo/AngularAnalysis/'];
% inputStruct(3).sample_prefix = 'U2OS_Xlone-FOXA2dDBD-Halo_';
% 
% inputStruct(4).dataPath = [dataPath_prefix '2021/20210403-PA646-U2OS-XLONE-FOXA2dCTD-Halo/AngularAnalysis/'];
% inputStruct(4).sample_prefix = 'U2OS_Xlone-FOXA2dCTD-Halo_';
% 
% inputStruct(5).dataPath = [dataPath_prefix '2021/20210429-PA646-U2OS-XLONE-FOXA2DBD-Halo/AngularAnalysis/'];
% inputStruct(5).sample_prefix = 'U2OS_Xlone-FOXA2DBD-Halo_';
% 
% inputStruct(6).dataPath = [dataPath_prefix '2022/20220412_U2OS_L30_OCT4-Halo/AngularAnalysis/'];
% inputStruct(6).sample_prefix = 'U2OS_L30_OCT4-Halo_';

for sampleIter = 1:length(inputStruct)
    frmRate = ["133Hz","74Hz"];
    for i = 1:length(frmRate)
        HMM_folder = fullfile(inputStruct(sampleIter).dataPath, 'HMM_first_QC_data');
        %find all MTT files:
        MTT_files=dir([HMM_folder, filesep, '*.mat']);
        temp_names = {MTT_files.name};
        % Filenames = ''; %for saving the actual file name
        idx = contains({MTT_files.name},frmRate(i));
        Filenames = temp_names{idx};
        HMM(i) = load(fullfile(MTT_files(1).folder,Filenames),'CellTracks','CellTrackViterbiClass');
    end

    % Merge the CellTracks
    CellTrackViterbiClass = horzcat(HMM(:).CellTrackViterbiClass);
    CellTracks = horzcat(HMM(:).CellTracks);
    
    savePath = fullfile(inputStruct(sampleIter).dataPath,'2TimePoints_mergeHMM_first_QC_data');
    if ~exist(savePath,'dir')
        mkdir(savePath)
    end
    
    save(fullfile(savePath,[inputStruct(sampleIter).sample_prefix  '133and74Hz_pooled_QC_CD2_classified.mat']),...
        'CellTracks', 'CellTrackViterbiClass');
end

%% Part II Local anisotropy calculation
close all; clear; clc;

dataPath_prefix = '/home/dell/Documents/METHOD/denglabpku-repository/anisotropy/ExampleData/';
FreSeg_size = [2 3 4 5 6]; % Required free segment length to do calculation
MaxJump = 0.85; % max trajector segment lengths allowed to calculate angles and anisotropy, unit: um
MinMinJumpThres = 0.100; % free segment smaller than this will not be analyzed, unit: um


inputStruct(1).dataPath = [dataPath_prefix 'AngularAnalysis/GlobalAnisotropy/2TimePoints_mergeHMM_first_QC_data/'];
inputStruct(1).dataName = 'U2OS_FOXA2-Halo_133and74Hz_pooled_QC_CD2_classified.mat';
inputStruct(1).sample_prefix = "U2OS FOXA2-Halo 133&74Hz";
inputStruct(1).MaxJump = MaxJump; 

% inputStruct(2).dataPath = [dataPath_prefix 'LocalAnisotropy\HMM_first_QC_data\'];
% inputStruct(2).dataName = 'U2OS_Xlone-FOXA2dNTD-Halo_133and74Hz_pooled_QC_CD2_classified.mat';
% inputStruct(2).sample_prefix = "U2OS ΔNTD-Halo 133&74Hz";
% inputStruct(2).MaxJump = MaxJump;
% 
% inputStruct(3).dataPath = [dataPath_prefix 'LocalAnisotropy\HMM_first_QC_data\'];
% inputStruct(3).dataName = 'U2OS_Xlone-FOXA2dDBD-Halo_133and74Hz_pooled_QC_CD2_classified.mat';
% inputStruct(3).sample_prefix = "U2OS ΔDBD-Halo 133&74Hz";
% inputStruct(3).MaxJump = MaxJump; 
% 
% inputStruct(4).dataPath = [dataPath_prefix 'LocalAnisotropy\HMM_first_QC_data\'];
% inputStruct(4).dataName = 'U2OS_Xlone-FOXA2dCTD-Halo_133and74Hz_pooled_QC_CD2_classified.mat';
% inputStruct(4).sample_prefix = "U2OS ΔCTD-Halo 133&74Hz";
% inputStruct(4).MaxJump = MaxJump; 
% 
% inputStruct(5).dataPath = [dataPath_prefix 'LocalAnisotropy\HMM_first_QC_data\'];
% inputStruct(5).dataName = 'U2OS_Xlone-FOXA2DBD-Halo_133and74Hz_pooled_QC_CD2_classified.mat';
% inputStruct(5).sample_prefix = "U2OS DBD-Halo 133&74Hz";
% inputStruct(5).MaxJump = MaxJump; 
% 
% inputStruct(6).dataPath = [dataPath_prefix 'LocalAnisotropy\HMM_first_QC_data\'];
% inputStruct(6).dataName = 'U2OS_L30_OCT4-Halo_133and74Hz_pooled_QC_CD2_classified.mat';
% inputStruct(6).sample_prefix = "U2OS OCT4-Halo 133&74Hz";
% inputStruct(6).MaxJump = MaxJump; 

% Start parallel free segments sorting and calculation
tic;

for i = 1:length(FreSeg_size)
    fprintf('=============== Current maximum free segment size: %d ==============\n',FreSeg_size(i));
    temp_FreSeg_size = FreSeg_size(i);
    parfor SampleIter = 1:length(inputStruct)
        inputStruct(SampleIter).BoundState = 1; % Designate bound state label, do not change
        inputStruct(SampleIter).FreeState = 2; % Designate free state label, do not change
        inputStruct(SampleIter).FreSeg_size = temp_FreSeg_size; % Required free segment length to do calculation
        inputStruct(SampleIter).nBins = 1*36; % number of bins to pool angles, do not change
        inputStruct(SampleIter).JackKnife_fraction = 0.5;  % JackKnife bootstrapping fraction, do not change
        inputStruct(SampleIter).JackKnife_iterations = 50; % JackKnife bootstrapping iteration, do not change
        inputStruct(SampleIter).MinMinJumpThres = MinMinJumpThres; 

        outputStruct = FreeSegmentAsymCalculation(inputStruct(SampleIter));

        FinalResults(i,SampleIter).samplePath = [inputStruct(SampleIter).dataPath inputStruct(SampleIter).dataName];
        FinalResults(i,SampleIter).sample_prefix = inputStruct(SampleIter).sample_prefix;
        FinalResults(i,SampleIter).FreSeg_size = temp_FreSeg_size;
        FinalResults(i,SampleIter).Alla = outputStruct.Alla; % theta(angle) for polarplot from all angles from pooled free segments after and before a bound segment
        FinalResults(i,SampleIter).AllVals = outputStruct.AllVals; % radius for polarplot from all angles from pooled free segments after and before a bound segment
        FinalResults(i,SampleIter).FreAfBo_mean_f_180_0 = outputStruct.FreAfBo_mean_f_180_0; % mean of f_180_0 of pooled free segments after a bound segment
        FinalResults(i,SampleIter).FreAfBo_std_f_180_0 = outputStruct.FreAfBo_std_f_180_0; % std of f_180_0 of pooled free segments after a bound segment
        FinalResults(i,SampleIter).FreBfBo_mean_f_180_0 = outputStruct.FreBfBo_mean_f_180_0; % mean of f_180_0 of pooled free segments before a bound segment
        FinalResults(i,SampleIter).FreBfBo_std_f_180_0 = outputStruct.FreBfBo_std_f_180_0; % std of f_180_0 of pooled free segments after a bound segment
        FinalResults(i,SampleIter).FreABfBo_mean_f_180_0 = outputStruct.FreABfBo_mean_f_180_0; % mean of f_180_0 of pooled free segments after and before a bound segment
        FinalResults(i,SampleIter).FreABfBo_std_f_180_0 = outputStruct.FreABfBo_std_f_180_0; % std of f_180_0 of pooled free segments after and before a bound segment
        FinalResults(i,SampleIter).FreeSeg_AfBo = outputStruct.FreeSeg_AfBo;% the index of Free segment after bound state segment, same as 'BoUnbo_CellTrackViterbiClass'
        FinalResults(i,SampleIter).FreeSeg_BfBo = outputStruct.FreeSeg_BfBo;% the index of Free segment before bound state segment
        FinalResults(i,SampleIter).BoSegAfterFreeSeg_BfB = outputStruct.BoSegAfterFreeSeg_BfB; % Consecutive bound segments after the free segments (those defined in FreeSeg_BfBo)
        FinalResults(i,SampleIter).IdxLUT_FreeSeg_AfBo = outputStruct.IdxLUT_FreeSeg_AfBo; % the same length and correspondingly to each element of 'FreeSeg_AfBo', but tell you the index of track index in 'BoUnbo_CellTrackViterbiClass'
        FinalResults(i,SampleIter).IdxLUT_FreeSeg_BfBo = outputStruct.IdxLUT_FreeSeg_BfBo; % the same length and correspondingly to each element of 'FreeSeg_BfBo', but tell you the index of track index in 'BoUnbo_CellTrackViterbiClass'
        FinalResults(i,SampleIter).angleMatrix_FreAfBo = outputStruct.angleMatrix_FreAfBo; % table contains free segments after bound: index of track index of BoUnboTrack (same as IdxLUT_FreeSeg_AfBo), first segment length, second segment length, angles in rad, complementary angles in rad
        FinalResults(i,SampleIter).angleMatrix_FreBfBo = outputStruct.angleMatrix_FreBfBo; % table contains free segments before bound: index of track index of BoUnboTrack (same as IdxLUT_FreeSeg_BfBo), first segment length, second segment length, angles in rad, complementary angles in rad
    end
end


fprintf('=== ANALYSIS COMPLETED! ===\n');
toc; 

% Save the result
result_output = 'LocalAnisotropy';
if ~exist([dataPath_prefix result_output],'dir')
    mkdir([dataPath_prefix result_output])
end
result_name = 'anisotropy_maxFreeAfBfbound_FOXA2related';
save([dataPath_prefix result_output filesep result_name '_' char(datetime('now','format','yyyy-MM-dd_HH_mm_ss')) '.mat']); 










