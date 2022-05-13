%   CompileTemporalSubSamplesOfHMM.m
%   written by Anders Sejr Hansen (AndersSejrHansen@post.harvard.edu;
%   @Anders_S_Hansen; https://anderssejrhansen.wordpress.com)
%   License: GNU GPL v3
%   Dependent functions:
%       - TemporallyReSampleCellTracks.m
%       - InferFrameRateFromName.m
clearvars -except ParentPath SampleNamePrefix resultDir; clc; close all;
% add path for dependent functions:
% addpath('Functions');

%   DESCRIPTION:
%   This script takes HMM-classified SPT data at high frame rates (223 Hz, 
%   133 Hz and 74 Hz) and then temporally sub-samples it to get data at
%   longer lag time. E.g. if you have data at 100 Hz and you take only the
%   data at odd frames (frame 1,3,5,7...) you end up getting data at 50 Hz.
%   One small concern when temporally subsampling is that when you
%   sub-sample to very long lag times (e.g. for data at 4.5 ms, 90 ms lag
%   time corresponds to frames 1,21,42,63,...) you end up with a slightly
%   biased set of trajectories, since the bound populations stays in focus
%   and gives long trajectories, whereas the free population largely
%   defocalizes. Thus, this could cause the HMM-classification to become
%   less robust and also means that slightly different trajectories could
%   be use for the different lag times.
%   To circumvent this, we do the HMM-classification first. First, the 223,
%   133 and 74 Hz data will be HMM-classified and then we will subsampled these
%   datasets. This is be more robust, although it will cost when it
%   comes to number of trajectories. The vbSPT-HMM cannot handle trajectory
%   gaps, so we will have to split into shorter trajectories, any
%   trajectory with gaps, which will cause us to have fewer very long
%   trajectories. But such is life. 

%   Modified by Zuhui 20200601 to analysis Wulan's fastSPT of APS, DE and H2B

% Script modified by Zuhui Wang 10/08/2020:
% Cell line: pl30-OCT4-DBD-FOXA2-CTD-detHelix
% Staining: Halo-PA646
% Acquisition: Nikon Ti2E, 7.5 ms, Left-camera, FF01-676/37.
% Comments: testing purpose, data collected on 09/31/2020, using home-made
% glass bottom 36-cm dish.

TempSubSampleStruc = struct();
SavePath = [resultDir, filesep, 'HMM_first_QC_data', filesep];

% Define the temporal subsamples: fps is calculated from mean frame
% duration.

% 15 ms (~66.7 Hz)
TempSubSampleStruc(1,1).TimePoints_conditions = {'223 Hz', '133 Hz', '74 Hz'};
TempSubSampleStruc(1,1).TimePointsToUse = cell(1,3);
TempSubSampleStruc(1,1).TimePointsToUse{1,1} = [];
TempSubSampleStruc(1,1).TimePointsToUse{1,2} = 2;
TempSubSampleStruc(1,1).TimePointsToUse{1,3} = [];
TempSubSampleStruc(1,1).SumSampleName = '66.7Hz';

% 18-22.5 ms (~44.4 Hz)
TempSubSampleStruc(1,2).TimePoints_conditions = {'223 Hz', '133 Hz', '74 Hz'};
TempSubSampleStruc(1,2).TimePointsToUse = cell(1,3);
TempSubSampleStruc(1,2).TimePointsToUse{1,1} = 4;
TempSubSampleStruc(1,2).TimePointsToUse{1,2} = 3;
TempSubSampleStruc(1,2).TimePointsToUse{1,3} = [];
TempSubSampleStruc(1,2).SumSampleName = '44.4Hz';

% 27-31.5 ms (~34.2 Hz)
TempSubSampleStruc(1,3).TimePoints_conditions = {'223 Hz', '133 Hz', '74 Hz'};
TempSubSampleStruc(1,3).TimePointsToUse = cell(1,3);
TempSubSampleStruc(1,3).TimePointsToUse{1,1} = [6,7];
TempSubSampleStruc(1,3).TimePointsToUse{1,2} = 4;
TempSubSampleStruc(1,3).TimePointsToUse{1,3} = 2;
TempSubSampleStruc(1,3).SumSampleName = '34.2Hz';

% 36-40.5 ms (~26 Hz)
TempSubSampleStruc(1,4).TimePoints_conditions = {'223 Hz', '133 Hz', '74 Hz'};
TempSubSampleStruc(1,4).TimePointsToUse = cell(1,3);
TempSubSampleStruc(1,4).TimePointsToUse{1,1} = [8,9];
TempSubSampleStruc(1,4).TimePointsToUse{1,2} = 5;
TempSubSampleStruc(1,4).TimePointsToUse{1,3} = 3;
TempSubSampleStruc(1,4).SumSampleName = '26Hz';

% 45 ms (~22.2 Hz)
TempSubSampleStruc(1,5).TimePoints_conditions = {'223 Hz', '133 Hz', '74 Hz'};
TempSubSampleStruc(1,5).TimePointsToUse = cell(1,3);
TempSubSampleStruc(1,5).TimePointsToUse{1,1} = 10;
TempSubSampleStruc(1,5).TimePointsToUse{1,2} = 6;
TempSubSampleStruc(1,5).TimePointsToUse{1,3} = [];
TempSubSampleStruc(1,5).SumSampleName = '22.2Hz';

% 52.5-54 ms (~18.8 Hz)
TempSubSampleStruc(1,6).TimePoints_conditions = {'223 Hz', '133 Hz', '74 Hz'};
TempSubSampleStruc(1,6).TimePointsToUse = cell(1,3);
TempSubSampleStruc(1,6).TimePointsToUse{1,1} = [11,12];
TempSubSampleStruc(1,6).TimePointsToUse{1,2} = 7;
TempSubSampleStruc(1,6).TimePointsToUse{1,3} = 4;
TempSubSampleStruc(1,6).SumSampleName = '18.8Hz';

% 58.5-63 ms (~16.5 Hz)
TempSubSampleStruc(1,7).TimePoints_conditions = {'223 Hz', '133 Hz', '74 Hz'};
TempSubSampleStruc(1,7).TimePointsToUse = cell(1,3);
TempSubSampleStruc(1,7).TimePointsToUse{1,1} = [13,14];
TempSubSampleStruc(1,7).TimePointsToUse{1,2} = 8;
TempSubSampleStruc(1,7).TimePointsToUse{1,3} = [];
TempSubSampleStruc(1,7).SumSampleName = '16.5Hz';

% 67.5 ms (~14.8 Hz)
TempSubSampleStruc(1,8).TimePoints_conditions = {'223 Hz', '133 Hz', '74 Hz'};
TempSubSampleStruc(1,8).TimePointsToUse = cell(1,3);
TempSubSampleStruc(1,8).TimePointsToUse{1,1} = 15;
TempSubSampleStruc(1,8).TimePointsToUse{1,2} = 9;
TempSubSampleStruc(1,8).TimePointsToUse{1,3} = 5;
TempSubSampleStruc(1,8).SumSampleName = '14.8Hz';

% % 81-82.5 ms (~12.2 Hz)
% TempSubSampleStruc(1,9).TimePoints_conditions = {'223 Hz', '133 Hz', '74 Hz'};
% TempSubSampleStruc(1,9).TimePointsToUse = cell(1,3);
% TempSubSampleStruc(1,9).TimePointsToUse{1,1} = [18,19];
% TempSubSampleStruc(1,9).TimePointsToUse{1,2} = 11;
% TempSubSampleStruc(1,9).TimePointsToUse{1,3} = 6;
% TempSubSampleStruc(1,9).SumSampleName = '12.2Hz';
% 
% % 90-97.5 ms (~10.6 Hz)
% TempSubSampleStruc(1,10).TimePoints_conditions = {'223 Hz', '133 Hz', '74 Hz'};
% TempSubSampleStruc(1,10).TimePointsToUse = cell(1,3);
% TempSubSampleStruc(1,10).TimePointsToUse{1,1} = [21,22];
% TempSubSampleStruc(1,10).TimePointsToUse{1,2} = [12,13];
% TempSubSampleStruc(1,10).TimePointsToUse{1,3} = 7;
% TempSubSampleStruc(1,10).SumSampleName = '10.6Hz';
% 
% % 105-112.5 ms (~9.2 Hz)
% TempSubSampleStruc(1,11).TimePoints_conditions = {'223 Hz', '133 Hz', '74 Hz'};
% TempSubSampleStruc(1,11).TimePointsToUse = cell(1,3);
% TempSubSampleStruc(1,11).TimePointsToUse{1,1} = [23,24,25];
% TempSubSampleStruc(1,11).TimePointsToUse{1,2} = [14,15];
% TempSubSampleStruc(1,11).TimePointsToUse{1,3} = 8;
% TempSubSampleStruc(1,11).SumSampleName = '9.2Hz';


% Compile sub-samples for all:
for DataSetIter = 1:1
    DataSet = DataSetIter;
    % Define the DataSets:
    if DataSet == 1
        MAT_files = struct();
        MAT_files.path = [resultDir, '/HMM_first_QC_data/'];
        MAT_files.Conditions = {'223 Hz', '133 Hz', '74 Hz'};
        MAT_files.SampleNamePrefix = SampleNamePrefix;
        MAT_files.workspaces = {'', ...
                                [MAT_files.SampleNamePrefix '133Hz_pooled_QC_CD2'],...
                                [MAT_files.SampleNamePrefix '74Hz_pooled_QC_CD2']};   
        MAT_files.SampleNamesuffix = '_pooled_QC_CD2_classified';
    end

    disp(['Processing DataSet number ', num2str(DataSetIter), ' of 1']);
disp('================================================');
    tic;
    % loop over all of the temporal subsamples to slice out:
    for iter = 1:length(TempSubSampleStruc)
        
        disp(['Proceesing temporal subsample ', num2str(iter), ' of ', num2str(length(TempSubSampleStruc)), ' total']);
        % number of movies sampled:
        MovieNumbCounter = 0;
        % make new placeholder cell arrays
        clear CellTracks CellTrackViterbiClass temp_CellTracks temp_CellTrackViterbiClass 
        temp_CellTracks = cell(1,1);
        temp_CellTrackViterbiClass = cell(1,1);
        
        %% This is a critical step to control which raw data set you use to subsampling
        % Now go through the 133 Hz and 74 Hz data one by one:
        for TimeIter = 2:length(TempSubSampleStruc(1,iter).TimePointsToUse)

            % check to see if there were any TimePoints to use:
            if ~isempty(TempSubSampleStruc(1,iter).TimePointsToUse{1,TimeIter})

                % Define temporal subsampling stuff:
                curr_TimePointsToUse_vector = TempSubSampleStruc(1,iter).TimePointsToUse{1,TimeIter};
                
                % Define workspace stuff
                curr_workspaces = MAT_files.workspaces{1, TimeIter};

                % load in each DataSet:           
                load([MAT_files.path, curr_workspaces, '_classified.mat']);
                
                % loop over each TimePointToUse in case there is more than one
                for TP_iter = 1:length(curr_TimePointsToUse_vector)
                    curr_TimePointsToUse = curr_TimePointsToUse_vector(TP_iter);

                    

                    % feed current dataset into the 
                    [ new_CellTracks, new_CellTrackViterbiClass ] = TemporallyReSampleCellTracks(CellTracks, CellTrackViterbiClass, curr_TimePointsToUse);
                    % now update the placeholders:
                    if MovieNumbCounter == 0
                        temp_CellTracks = new_CellTracks;
                        temp_CellTrackViterbiClass = new_CellTrackViterbiClass;
                        MovieNumbCounter = MovieNumbCounter + 1;
                    else
                        % append the new trackedPar (append placeholder cell arrays):
                        temp_CellTracks = horzcat(temp_CellTracks, new_CellTracks);
                        temp_CellTrackViterbiClass = horzcat(temp_CellTrackViterbiClass, new_CellTrackViterbiClass);
                        MovieNumbCounter = MovieNumbCounter + 1;
                    end
                end
            end

        end

        clear CellTracks CellTrackViterbiClass
        CellTracks = temp_CellTracks;
        CellTrackViterbiClass = temp_CellTrackViterbiClass;
        
        % get the LagTime:
        LagTime = InferFrameRateFromName( [MAT_files.SampleNamePrefix, TempSubSampleStruc(1,iter).SumSampleName,  MAT_files.SampleNamesuffix, '.mat']);

        %%%% FINISHED A FULL COLLECTION FOR A NEW FRAME-RATE: SAVE THE DATA %%%
        %   just save all the data to a single trackedPar to avoid having a lot
        %   of files:
        disp(['saving file: ', [MAT_files.SampleNamePrefix, TempSubSampleStruc(1,iter).SumSampleName,  MAT_files.SampleNamesuffix], '...']);
        save([SavePath, MAT_files.SampleNamePrefix, TempSubSampleStruc(1,iter).SumSampleName,  MAT_files.SampleNamesuffix, '.mat'], 'CellTracks', 'CellTrackViterbiClass');
        
    end
    toc;
    disp('================================================');
end



