%   MergeQC_SPT_data.m
%   written by Anders Sejr Hansen (AndersSejrHansen@post.harvard.edu;
%   @Anders_S_Hansen; https://anderssejrhansen.wordpress.com)
%   License: GNU GPL v3
%   Dependent functions:
%       - RemoveAmbigiousTracks.m

clear; clc; close all;
%   DESCRIPTION
%   Merge and QC SPT data: 
%       - read in individual MAT files from single cells
%       - Filter out trajectories that got too close
%       - merge all individual trackedPar into a single long one

%% Modify the path but keep the rest
%% Before the start, remember to add the parent folder of anisotropy to matlab path
ParentPath = '/home/dell/Documents/METHOD/denglabpku-repository/anisotropy/ExampleData';
SampleNamePrefix = 'U2OS_FOXA2-Halo_';

%% PROCESS START
ClosestDist = 2; % closest distance between particle in micrometers. This is the user-specified threshold that determines which particles will be removed

inputPath = [ParentPath filesep 'MTT_trajectory']; % without filesep

% Make output folder
mkdir([ParentPath filesep 'AngularAnalysis' filesep 'GlobalAnisotropy' filesep 'QC_data']);
mkdir([ParentPath filesep 'AngularAnalysis' filesep 'GlobalAnisotropy' filesep 'QC_data_reformatted']);
mkdir([ParentPath filesep 'AngularAnalysis' filesep 'GlobalAnisotropy' filesep 'HMM_first_QC_data']);
mkdir([ParentPath filesep 'AngularAnalysis' filesep 'GlobalAnisotropy' filesep 'FolderWithRunInputFiles']);
mkdir([ParentPath filesep 'AngularAnalysis' filesep 'GlobalAnisotropy' filesep 'SpatioTemporalAnalysis']);
mkdir([ParentPath filesep 'AngularAnalysis' filesep 'GlobalAnisotropy' filesep 'QC_Plots']);

% where to save to?
resultDir = [ParentPath filesep 'AngularAnalysis' filesep 'GlobalAnisotropy' filesep];

%Define frame rates string contained in file name
frameRate = {'133Hz','74Hz'};
% SampleName

SampleName = {...
    [SampleNamePrefix '133Hz_pooled'],...
    [SampleNamePrefix '74Hz_pooled']};

%find all MTT files:
MTT_files=dir([inputPath, filesep,'*.mat']);
Filenames = {MTT_files.name}; %for saving the actual file name

for FrmRate_iter = 1:length(frameRate) %
    input_struct(FrmRate_iter).path = [inputPath, filesep];
    Idx = contains(Filenames, frameRate(FrmRate_iter));
    input_struct(FrmRate_iter).workspaces = Filenames(Idx);
end
    

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% PERFORM QC ON DATA %%%%%%%%%%%%%%%%%%%%%%%%
for FrmRate_iter = 1:length(input_struct)
    disp('==============================================================');
    tic;
    disp(['Current Analyzing Frame Rate: ', frameRate{FrmRate_iter}]);

    trackedPar_QC = struct();
    trackedParCounter = 0;
    % loop over each replicate

    curr_path = input_struct(FrmRate_iter).path;
    curr_workspaces = input_struct(FrmRate_iter).workspaces;

    for CellIter = 1:length(curr_workspaces)

            % load in the data
            load([curr_path, curr_workspaces{CellIter}]);
            % perform QC
            temp_trackedPar_merged = RemoveAmbigiousTracks( trackedPar, ClosestDist );

            % because you are silly, your old data uses a row vector
            % format for Frame and TimeStamp and your new format uses a
            % column vector for this; so need to check if they were in
            % a row vector format and then convert to column vector
            % format:
            for TrackIter = 1:length(temp_trackedPar_merged)
                if isrow(temp_trackedPar_merged(1,TrackIter).Frame)
                    temp_trackedPar_merged(1,TrackIter).Frame = temp_trackedPar_merged(1,TrackIter).Frame';
                    temp_trackedPar_merged(1,TrackIter).TimeStamp = temp_trackedPar_merged(1,TrackIter).TimeStamp';
                end
            end



            % merge the dataset
            if trackedParCounter == 0
                trackedPar_QC = temp_trackedPar_merged;
                % update the counter:
                trackedParCounter = trackedParCounter + 1;
            else
                % how much to add to the frame number
                ToAdd = trackedParCounter * settings.Frames; % * NumbFrames
                % update the frame number:
                for TrackIter = 1:length(temp_trackedPar_merged)
                    temp_trackedPar_merged(1,TrackIter).Frame = temp_trackedPar_merged(1,TrackIter).Frame + ToAdd;
                end
                % append the new trackedPar:
                trackedPar_QC = horzcat(trackedPar_QC, temp_trackedPar_merged);
                % update the counter:
                trackedParCounter = trackedParCounter + 1;
            end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%% SAVE THE QC'ed DATA %%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('saving the merged and QCed dataset %s \n',SampleName{FrmRate_iter});
    save([resultDir,'QC_data', filesep, SampleName{FrmRate_iter}, '_QC_CD', num2str(ClosestDist), '.mat'], 'trackedPar_QC');
%     end
    
    
    
    toc;
    disp('==============================================================');
end

fprintf('Starting run step 2...\n');
run Step2_Batch_vbSPT_classify_merge2time.m
fprintf('Starting run step 3...\n');
run Step3_CompileTemporalSubSamplesOfHMM_merge2time.m
fprintf('Starting run step 4...\n');
run Step4_Process_SpatioTemporal_AngleAnalysis_v2_merge2time.m
fprintf('Step 4 finished...\n');