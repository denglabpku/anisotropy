% Merge multiple SptioTemporalAnalysis results from different cell line
% into one FinalResults
clear; close all; clc;

saveDir = '/mnt/f58069a5-1cf3-43b8-bb9b-ea74327327c9/WZH-DataCenter/PROCESS-SPT/AnisotropyPlot_XLONEFOXA2related/';
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end


% XLONE-FOXA2-FL
FinalResult_path{1} = '/mnt/f58069a5-1cf3-43b8-bb9b-ea74327327c9/WZH-DataCenter/PROCESS-SPT/2021/20210424-PA646-U2OS-XLONE-FOXA2(D6)-Halo/AngularAnalysis/2TimePoints/SpatioTemporalAnalysis/';
FinalResult_name{1} = 'PA646_SpatioTemporalAnalysis_2022-04-14_11_44_35.mat';

% XLONE-FOXA2dNTD
FinalResult_path{2} = '/mnt/f58069a5-1cf3-43b8-bb9b-ea74327327c9/WZH-DataCenter/PROCESS-SPT/2021/20210413-PA646-U2OS-XLONE-FOXA2dNTD-Halo/AngularAnalysis/2TimePoints/SpatioTemporalAnalysis/';
FinalResult_name{2} = 'PA646_SpatioTemporalAnalysis_2022-04-14_14_39_06.mat';

% XLONE-FOXA2dDBD
FinalResult_path{3} = '/mnt/f58069a5-1cf3-43b8-bb9b-ea74327327c9/WZH-DataCenter/PROCESS-SPT/2021/20210506-PA646-U2OS-XLONE-FOXA2dDBD-Halo/AngularAnalysis/2TimePoints/SpatioTemporalAnalysis/';
FinalResult_name{3} = 'PA646_SpatioTemporalAnalysis_2022-04-14_14_28_09.mat';

% XLONE-FOXA2dCTD
FinalResult_path{4} = '/mnt/f58069a5-1cf3-43b8-bb9b-ea74327327c9/WZH-DataCenter/PROCESS-SPT/2021/20210403-PA646-U2OS-XLONE-FOXA2dCTD-Halo/AngularAnalysis/2TimePoints/SpatioTemporalAnalysis/';
FinalResult_name{4} = 'PA646_SpatioTemporalAnalysis_2022-04-14_15_19_26.mat';

% XLONE-FOXA2DBD
FinalResult_path{5} = '/mnt/f58069a5-1cf3-43b8-bb9b-ea74327327c9/WZH-DataCenter/PROCESS-SPT/2021/20210429-PA646-U2OS-XLONE-FOXA2DBD-Halo/AngularAnalysis/2TimePoints/SpatioTemporalAnalysis/';
FinalResult_name{5} = 'PA646_SpatioTemporalAnalysis_2022-04-14_14_33_00.mat';


for sampleIter = 1:length(FinalResult_name)
    load([FinalResult_path{sampleIter} FinalResult_name{sampleIter}],'FinalResults', 'input_struct');
    Merge_FinalResults(sampleIter) = FinalResults;
    Merge_input_struct(sampleIter) = input_struct;
    clear FinalResults input_struct
end

save([saveDir 'Merged2TimePoints_PA646_SpatioTemporalAnalysis_' char(datetime('now','format','yyyy-MM-dd_HH_mm_ss'))])

fprintf('Merging FinalResults completed! \n');
    