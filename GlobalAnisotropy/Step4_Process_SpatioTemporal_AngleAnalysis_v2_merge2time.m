
%   Process_SpatioTemporal_AngleAnalysis_v2.m
%   written by Anders Sejr Hansen (AndersSejrHansen@post.harvard.edu;
%   @Anders_S_Hansen; https://anderssejrhansen.wordpress.com)
%   License: GNU GPL v3
clearvars -except ParentPath SampleNamePrefix resultDir; clc; close all;

%   DESCRIPTION
%   Perform angle analysis for fastSPT data at multiple spatial and
%   temporal scales. In particular, for a given cell line/condition,
%   calculate and plot the angle anisotropy for across space (as a function
%   of jump length) and across time (as a function of FrameRate). 

%   UPDATE Version 2 - October 26, 2017
%   Updated to add a few new features:
%   -   count number of angles:
%       - for all jumps
%       - per displacement length bin
%   This will be helpful for taking weighted averages
%   - set global min-min threshold: only if jumps exceed this threshold
%   will they be counted. This will be helpful for getting rid of very
%   small angles. 
%   Do weighted average

% Script modified by Zuhui Wang 10/08/2020:
% Cell line: pl30-OCT4-DBD-FOXA2-CTD-detHelix
% Staining: Halo-PA646
% Acquisition: Nikon Ti2E, 7.5 ms, Left-camera, FF01-676/37.
% Comments: testing purpose, data collected on 09/31/2020, using home-made
% glass bottom 36-cm dish.

% add path for dependent functions:
% addpath('Functions');
saveDir = fullfile(resultDir, 'SpatioTemporalAnalysis');
data_path = [resultDir filesep 'HMM_first_QC_data' filesep];
DataSet = 1;
datetime = char(datetime('now','format','yyyy-MM-dd_HH_mm_ss'));

% Define which DataSet to analyze:
sample_suffix = '_pooled_QC_CD2_classified.mat';

if DataSet == 1
    sample_prefixes = {SampleNamePrefix}; % add more if analyzing more than one dataset
    SaveMatName = 'PA646_SpatioTemporalAnalysis_'; % name of the final file to be saved
end


% 222 Hz
% Analysis_struct(1,1).FrameRate = '223Hz';
% Analysis_struct(1,1).LagTime = 0.0045;
% Analysis_struct(1,1).MaxJump = 0.45; %in μm
% 133 Hz
Analysis_struct(1,1).FrameRate = '133Hz';
Analysis_struct(1,1).LagTime = 0.0075;
Analysis_struct(1,1).MaxJump = 0.55;
% 74 Hz
Analysis_struct(1,2).FrameRate = '74Hz';
Analysis_struct(1,2).LagTime = 0.0135;
Analysis_struct(1,2).MaxJump = 0.65;
% 66.7 Hz
Analysis_struct(1,3).FrameRate = '66.7Hz';
Analysis_struct(1,3).LagTime = 0.015;
Analysis_struct(1,3).MaxJump = 0.70;
% 44.4 Hz
Analysis_struct(1,4).FrameRate = '44.4Hz';
Analysis_struct(1,4).LagTime = 1/44.4;
Analysis_struct(1,4).MaxJump = 0.75;
% 34.2 Hz
Analysis_struct(1,5).FrameRate = '34.2Hz';
Analysis_struct(1,5).LagTime = 1/34.2;
Analysis_struct(1,5).MaxJump = 1.05;
% 26 Hz
Analysis_struct(1,6).FrameRate = '26Hz';
Analysis_struct(1,6).LagTime = 1/26;
Analysis_struct(1,6).MaxJump = 1.2;
% 22.2 Hz
Analysis_struct(1,7).FrameRate = '22.2Hz';
Analysis_struct(1,7).LagTime = 1/22.2;
Analysis_struct(1,7).MaxJump = 1.2;
% 18.8 Hz
Analysis_struct(1,8).FrameRate = '18.8Hz';
Analysis_struct(1,8).LagTime = 1/18.8;
Analysis_struct(1,8).MaxJump = 1.2;
% 16.5 Hz
Analysis_struct(1,9).FrameRate = '16.5Hz';
Analysis_struct(1,9).LagTime = 1/16.5;
Analysis_struct(1,9).MaxJump = 1.2;
% 14.8 Hz
Analysis_struct(1,10).FrameRate = '14.8Hz';
Analysis_struct(1,10).LagTime = 1/14.8;
Analysis_struct(1,10).MaxJump = 1.2;
% % 12.2 Hz
% Analysis_struct(1,10).FrameRate = '12.2Hz';
% Analysis_struct(1,10).LagTime = 1/12.2;
% Analysis_struct(1,10).MaxJump = 1.2;
% % 10.6 Hz
% Analysis_struct(1,11).FrameRate = '10.6Hz';
% Analysis_struct(1,11).LagTime = 1/10.6;
% Analysis_struct(1,11).MaxJump = 1.2;
% % 9.2 Hz
% Analysis_struct(1,12).FrameRate = '9.2Hz';
% Analysis_struct(1,12).LagTime = 1/9.2;
% Analysis_struct(1,12).MaxJump = 1.2;

% Define structure for passing to the Angle Analysis function:
JackKnife_fraction = 0.5;
input_struct(1).JackKnife_fraction = JackKnife_fraction;
JackKnife_iterations = 50;
input_struct(1).JackKnife_iterations = JackKnife_iterations;
%%%%%%%%%%% KEY PARAMETERS %%%%%%%%%%%
%How many bins to use?
nBins = 1*36; %so 10 degrees per bin
input_struct(1).nBins = nBins;
rads = 0:2*pi/nBins:2*pi;
input_struct(1).rads = rads;
%Classify the states
BoundState = 1;
input_struct(1).BoundState = BoundState;
FreeState = 2;
input_struct(1).FreeState = FreeState;
minFree = 2;
minAngleNumber = 200;
input_struct(1).minAngleNumber = minAngleNumber; % ensure there are enough angle data, below this number, don't do calculations
input_struct(1).minFree = minFree; % at least 2 frames are required to compute an angle
MovingThreshold = 0.10:0.05:0.950; % range of translocations
MinMinJumpThres = 0.100; %This threshold is the minimum for both jumps making up the angle
input_struct(1).MinMinJumpThres = MinMinJumpThres;

% Set a GLOBAL threshold for minimal displacements that will apply even to
% analysis at multiple spatiotemporal scales:
GlobalMinJumpThres = 0.1; % so 100 nm
input_struct(1).GlobalMinJumpThres = GlobalMinJumpThres;

% define a threshold to reduce strang long trajectory artifact
MinNumAngles = 5;
MaxAsymAnglesFrac = 0.5;
input_struct(1).MinNumAngles = MinNumAngles;
input_struct(1).MaxAsymAnglesFrac = MaxAsymAnglesFrac;

% threshold for calculating offset: 0-72 degrees; 0-1.26 rads
% See SupMaterial of CTCF NatChemBio for more detail
% Briefly, 72 degress is considered to be background of isotropy and
% substracted from the angle PDF and to see the anisotropy effect.
OffSet_thres = 1.26;
input_struct(1).OffSet_thres = OffSet_thres;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% PERFORM ANALYSIS %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TotalDataToAnalyze = length(sample_prefixes) * length(Analysis_struct);
DataIter = 1;
% loop over the samples
for SampleIter = 1:length(sample_prefixes)
    disp('==========================================================');
    disp(['Analyzing ',  sample_prefixes{SampleIter}]);
    disp(['Analyzing sample number ', num2str(SampleIter), ' of ', num2str(length(sample_prefixes)), ' total samples']);
    % assign a name to the current sample:
    FinalResults(1,SampleIter).SampleName = sample_prefixes{SampleIter};
    % initialize vectors, matrices and structures to handle the output
    % data:
    FinalResults(1,SampleIter).FrameRates = Analysis_struct;
    % for saving distributions:
    FinalResults(1,SampleIter).Alla = cell(1,length(Analysis_struct));
    FinalResults(1,SampleIter).AllVals = cell(1,length(Analysis_struct));
    FinalResults(1,SampleIter).normPDF = cell(1,length(Analysis_struct));
    % for saving overall metrics:
    FinalResults(1,SampleIter).mean_Amp = zeros(1,length(Analysis_struct));
    FinalResults(1,SampleIter).mean_FWHM = zeros(1,length(Analysis_struct));
    FinalResults(1,SampleIter).mean_AC = zeros(1,length(Analysis_struct));
    FinalResults(1,SampleIter).mean_f_180_0 = zeros(1,length(Analysis_struct));
    FinalResults(1,SampleIter).TotalNumAngles = zeros(1,length(Analysis_struct));
    FinalResults(1,SampleIter).TotalMinMinAngles = zeros(1,length(Analysis_struct));
    % for saving Jack-Knife sampled data
    FinalResults(1,SampleIter).jack_Amp = cell(1,length(Analysis_struct));
    FinalResults(1,SampleIter).jack_FWHM = cell(1,length(Analysis_struct));
    FinalResults(1,SampleIter).jack_AC = cell(1,length(Analysis_struct));
    FinalResults(1,SampleIter).jack_f_180_0 = cell(1,length(Analysis_struct));
    % for saving Spatio-Temporal Matrices (overall means):
        % 1st dimension: Frame Rate (time)
        % 2nd dimension: displacement length (space)
        FinalResults(1,SampleIter).AmpAtClosestMean = NaN(length(Analysis_struct), length(MovingThreshold));
        FinalResults(1,SampleIter).FWHMAtClosestMean = NaN(length(Analysis_struct), length(MovingThreshold));
        FinalResults(1,SampleIter).ACAtClosestMean = NaN(length(Analysis_struct), length(MovingThreshold));
        FinalResults(1,SampleIter).f_180_0_AtClosestMean = NaN(length(Analysis_struct), length(MovingThreshold));
        FinalResults(1,SampleIter).AmpAtClosestMin = NaN(length(Analysis_struct), length(MovingThreshold));
        FinalResults(1,SampleIter).FWHMAtClosestMin = NaN(length(Analysis_struct), length(MovingThreshold));
        FinalResults(1,SampleIter).ACAtClosestMin = NaN(length(Analysis_struct), length(MovingThreshold));
        FinalResults(1,SampleIter).f_180_0_AtClosestMin = NaN(length(Analysis_struct), length(MovingThreshold));
        FinalResults(1,SampleIter).NumAnglesAtClosestMean = NaN(length(Analysis_struct), length(MovingThreshold));
        FinalResults(1,SampleIter).NumAnglesAtClosestMin = NaN(length(Analysis_struct), length(MovingThreshold));
    % for saving Spatio-Temporal Matrices (overall means):
        % use a cell array of vectors (i.e. a matrix of vectors):
        % 1st dimension: Frame Rate (time)
        % 2nd dimension: displacement length (space)
        % 3rd dimension: different iterations of Jack-Knife sampling:
        FinalResults(1,SampleIter).jack_AmpAtClosestMean = cell(length(Analysis_struct), length(MovingThreshold));
        FinalResults(1,SampleIter).jack_FWHMAtClosestMean = cell(length(Analysis_struct), length(MovingThreshold));
        FinalResults(1,SampleIter).jack_ACAtClosestMean = cell(length(Analysis_struct), length(MovingThreshold));
        FinalResults(1,SampleIter).jack_f_180_0_AtClosestMean = cell(length(Analysis_struct), length(MovingThreshold));
        FinalResults(1,SampleIter).jack_AmpAtClosestMin = cell(length(Analysis_struct), length(MovingThreshold));
        FinalResults(1,SampleIter).jack_FWHMAtClosestMin = cell(length(Analysis_struct), length(MovingThreshold));
        FinalResults(1,SampleIter).jack_ACAtClosestMin = cell(length(Analysis_struct), length(MovingThreshold));
        FinalResults(1,SampleIter).jack_f_180_0_AtClosestMin = cell(length(Analysis_struct), length(MovingThreshold));
    % for saving the full matrix of angles vs. jump lengths:
    FinalResults(1,SampleIter).CellOfAngleMatrix = cell(1, length(Analysis_struct));
    
    % Loop over the frame rates
    for TimeIter = 1:length(Analysis_struct)
        tic;
        disp(['    Analyzing frame rate ', num2str(TimeIter), ' of ', num2str(length(Analysis_struct)), ' total']);
        disp(['    Analyzing dataset #', num2str(DataIter), ' of ', num2str(TotalDataToAnalyze), ' datasets']);
        
        %%%%%%%% STEP 1: Load in the current data %%%%%%%%
        clear CellTracks CellTrackViterbiClass
        load([data_path, sample_prefixes{SampleIter}, Analysis_struct(1,TimeIter).FrameRate, sample_suffix]);
        % add the relevant data
        input_struct(1).CellTracks = CellTracks;
        input_struct(1).CellTrackViterbiClass = CellTrackViterbiClass;
        % set an upper limit on how big a jump will be considered:
        input_struct(1).MaxJump = Analysis_struct(1,TimeIter).MaxJump + 0.3;
        % update the Max Jump to consider (for some frame rates, we want to limit how long displacements we consider):
        MaxJump = Analysis_struct(1,TimeIter).MaxJump;
        [~, CurrIndex] = min(abs(MovingThreshold - MaxJump)); 
        curr_MovingThreshold = MovingThreshold(1:(CurrIndex-1));
        input_struct(1).MovingThreshold = curr_MovingThreshold; % update Moving threshold
        
        %%%%%%%% STEP 2: PERFORM ANGLE ANALYSIS ON CURRENT DATA %%%%%%%%
        % perform analysis for the current sample and frame rate:
        output_struct = angleFWHM_Amp_HMM_analyzer_v5(input_struct);
        
        
        %%%%%%% STEP 3: SAVE THE CURRENT DATA IN A (SENSIBLE) WAY  %%%%%%%
        FinalResults(1,SampleIter).FrameRates = Analysis_struct;
        % for saving distributions:
        FinalResults(1,SampleIter).Alla{1,TimeIter} = output_struct.Alla;
        FinalResults(1,SampleIter).AllVals{1,TimeIter} = output_struct.AllVals;
        FinalResults(1,SampleIter).normPDF{1,TimeIter} = output_struct.normPDF;
        % for saving overall metrics:
        FinalResults(1,SampleIter).mean_Amp(1,TimeIter) = output_struct.mean_Amp;
        FinalResults(1,SampleIter).mean_FWHM(1,TimeIter) = output_struct.mean_FWHM;
        FinalResults(1,SampleIter).mean_AC(1,TimeIter) = output_struct.AC;
        FinalResults(1,SampleIter).mean_f_180_0(1,TimeIter) = output_struct.f_180_0;
        FinalResults(1,SampleIter).TotalNumAngles(1,TimeIter) = output_struct.TotalNumAngles;
        FinalResults(1,SampleIter).TotalMinMinAngles(1,TimeIter) = output_struct.TotalMinMinAngles;
        % for saving Jack-Knife sampled data
        FinalResults(1,SampleIter).jack_Amp{1,TimeIter} = output_struct.jack_Amp;
        FinalResults(1,SampleIter).jack_FWHM{1,TimeIter} = output_struct.jack_FWHM;
        FinalResults(1,SampleIter).jack_AC{1,TimeIter} = output_struct.jack_AC;
        FinalResults(1,SampleIter).jack_f_180_0{1,TimeIter} = output_struct.jack_f_180_0;
        % save cell array of matrics:
        FinalResults(1,SampleIter).CellOfAngleMatrix{1,TimeIter} = output_struct.AngleMatrix;
        % for saving Spatio-Temporal Matrices (overall means):
            % 1st dimension: Frame Rate (time)
            % 2nd dimension: displacement length (space)
            FinalResults(1,SampleIter).AmpAtClosestMean(TimeIter, 1:length(output_struct.AmpAtClosestMean)) = output_struct.AmpAtClosestMean;
            FinalResults(1,SampleIter).FWHMAtClosestMean(TimeIter, 1:length(output_struct.FWHMAtClosestMean)) = output_struct.FWHMAtClosestMean;
            FinalResults(1,SampleIter).ACAtClosestMean(TimeIter, 1:length(output_struct.ACAtClosestMean)) = output_struct.ACAtClosestMean;
            FinalResults(1,SampleIter).f_180_0_AtClosestMean(TimeIter, 1:length(output_struct.f_180_0_AtClosestMean)) = output_struct.f_180_0_AtClosestMean;
            FinalResults(1,SampleIter).AmpAtClosestMin(TimeIter, 1:length(output_struct.AmpAtClosestMin)) = output_struct.AmpAtClosestMin;
            FinalResults(1,SampleIter).FWHMAtClosestMin(TimeIter, 1:length(output_struct.FWHMAtClosestMin)) = output_struct.FWHMAtClosestMin;
            FinalResults(1,SampleIter).ACAtClosestMin(TimeIter, 1:length(output_struct.ACAtClosestMin)) = output_struct.ACAtClosestMin;
            FinalResults(1,SampleIter).f_180_0_AtClosestMin(TimeIter, 1:length(output_struct.f_180_0_AtClosestMin)) = output_struct.f_180_0_AtClosestMin;
            % save the numbers:
            FinalResults(1,SampleIter).NumAnglesAtClosestMean(TimeIter, 1:length(output_struct.NumAnglesAtClosestMean)) = output_struct.NumAnglesAtClosestMean;
            FinalResults(1,SampleIter).NumAnglesAtClosestMin(TimeIter, 1:length(output_struct.NumAnglesAtClosestMin)) = output_struct.NumAnglesAtClosestMin;
        % for saving Spatio-Temporal Matrices (overall means):
            % use a cell array of vectors (i.e. a matrix of vectors):
            % 1st dimension: Frame Rate (time)
            % 2nd dimension: displacement length (space)
            % 3rd dimension: different iterations of Jack-Knife sampling:
            
            % loop over each bin in the moving threshold vector and save
            % the Jack-Knife sampling
            for CurrIter = 1:length(curr_MovingThreshold)
                FinalResults(1,SampleIter).jack_AmpAtClosestMean{TimeIter,CurrIter} = output_struct.jack_AmpAtClosestMean(:,CurrIter);
                FinalResults(1,SampleIter).jack_FWHMAtClosestMean{TimeIter,CurrIter} = output_struct.jack_FWHMAtClosestMean(:,CurrIter);
                FinalResults(1,SampleIter).jack_ACAtClosestMean{TimeIter,CurrIter} = output_struct.jack_ACAtClosestMean(:,CurrIter);
                FinalResults(1,SampleIter).jack_f_180_0_AtClosestMean{TimeIter,CurrIter} = output_struct.jack_f_180_0_AtClosestMean(:,CurrIter);
                FinalResults(1,SampleIter).jack_AmpAtClosestMin{TimeIter,CurrIter} = output_struct.jack_AmpAtClosestMin(:,CurrIter);
                FinalResults(1,SampleIter).jack_FWHMAtClosestMin{TimeIter,CurrIter} = output_struct.jack_FWHMAtClosestMin(:,CurrIter);
                FinalResults(1,SampleIter).jack_ACAtClosestMin{TimeIter,CurrIter} = output_struct.jack_ACAtClosestMin(:,CurrIter);
                FinalResults(1,SampleIter).jack_f_180_0_AtClosestMin{TimeIter,CurrIter} = output_struct.jack_f_180_0_AtClosestMin(:,CurrIter);
            end
         DataIter = DataIter + 1;
        clear output_struct
        toc;
    end
    
    disp('calculating the full matrix across all displacement length bins'); tic;
    % For the full AngleMatrix as a function of displacement length,
    % average over all timepoints:
    [ Amp_JumpMatrix, FWHM_JumpMatrix, AC_JumpMatrix, f_180_0_JumpMatrix,  NumAngles_JumpMatrix ] = AngleMatrix_analyzer( input_struct, FinalResults(1,SampleIter).CellOfAngleMatrix );
    % save the key output
    FinalResults(1,SampleIter).Amp_JumpMatrix = Amp_JumpMatrix;
    FinalResults(1,SampleIter).FWHM_JumpMatrix = FWHM_JumpMatrix;
    FinalResults(1,SampleIter).AC_JumpMatrix = AC_JumpMatrix;
    FinalResults(1,SampleIter).f_180_0_JumpMatrix = f_180_0_JumpMatrix;
    FinalResults(1,SampleIter).NumAngles_JumpMatrix = NumAngles_JumpMatrix;
    toc;
    
    
    disp('==========================================================');
    
end
disp('=================Export to mat file========================');
save([saveDir, filesep, SaveMatName, datetime, '.mat'], 'FinalResults', 'input_struct','datetime');
disp('=================Export to mat file done========================');















