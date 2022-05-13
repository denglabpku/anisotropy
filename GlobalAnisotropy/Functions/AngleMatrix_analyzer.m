function [ Amp_JumpMatrix, FWHM_JumpMatrix, AC_JumpMatrix, f_180_0_JumpMatrix, NumAngles_JumpMatrix ] = AngleMatrix_analyzer( input_struct, CellOfAngleMatrix )
%ANGLEMATRIX_ANALYZER Summary of this function goes here
%   Detailed explanation goes here
% conveniently store the inputs
nBins = input_struct(1).nBins;
rads = input_struct(1).rads;
MinNumAngles = input_struct(1).MinNumAngles;
OffSet_thres = input_struct(1).OffSet_thres;


% perform the analysis in two steps:
%   1st: loop over all timepoints and make a single unified cell array of
%   matrix of angles
%   2nd: loop over each displacement bin and calculate Amp, AC, FWHM,
%   f(180/0)

%%%%%%%%%%% STEP 1: make unified matrix: %%%%%%%%%%%
max_dim = 0;
for iter=1:length(CellOfAngleMatrix)
    if max(size(CellOfAngleMatrix{1,iter})) > max_dim
        max_dim = max(size(CellOfAngleMatrix{1,iter}));
    end
end
MergedAngleMatrix = cell(max_dim,max_dim);
% now merge, but be careful because the sub-cell arrays have different
% sizes:
for iter=1:length(CellOfAngleMatrix)
    curr_cell = CellOfAngleMatrix{1,iter};
    for RowIter = 1:size(curr_cell,1)
        for ColIter = 1:size(curr_cell,2)
            MergedAngleMatrix{RowIter, ColIter} = horzcat(MergedAngleMatrix{RowIter, ColIter}, curr_cell{RowIter, ColIter});
        end
    end
end

%%%%%%%%%%% STEP 2: Perform angle calculations %%%%%%%%%%%
% calculate Amplitude, AC, FWHM, f_180_0 for the full matrix:
Amp_JumpMatrix = zeros(max_dim, max_dim);
FWHM_JumpMatrix = zeros(max_dim, max_dim);
AC_JumpMatrix = zeros(max_dim, max_dim);
f_180_0_JumpMatrix = zeros(max_dim, max_dim);
NumAngles_JumpMatrix = zeros(max_dim, max_dim);

% calculate Amplitude, AC, FWHM, f_180_0 for the full matrix:
for RowIter = 1:max_dim
    for ColIter = 1:max_dim
        %%%%%%%%%%%% CALCULATION FOR FULL MATRIX %%%%%%%%%%%%
        % ensure that there is enough data
        if numel(MergedAngleMatrix{RowIter, ColIter}) > MinNumAngles
            [Alla,AllVals] = rose(MergedAngleMatrix{RowIter, ColIter}, nBins);
            AllVals = AllVals./sum(AllVals);
            %Calculate the Assymmetry Coefficient (Izeddin et al.)
            %Use 0 +/- 30 degrees and 180 +/- 30 degrees as the 
            NonRedundantVals = AllVals(2:4:end);
            NumElements = round(length(NonRedundantVals)*30/360);
            For = [1:1:NumElements length(NonRedundantVals)-NumElements:1:length(NonRedundantVals)];
            Back = [round(length(NonRedundantVals)/2)-NumElements+1:1:round(length(NonRedundantVals)/2)+NumElements];
            AC_JumpMatrix(RowIter, ColIter) = log2(mean(NonRedundantVals(For))/mean(NonRedundantVals(Back)));
            f_180_0_JumpMatrix(RowIter, ColIter) = mean(NonRedundantVals(Back))/mean(NonRedundantVals(For));
            % make a normalised PDF:
            normPDF = NonRedundantVals ./sum(NonRedundantVals);
            [ Amp_JumpMatrix(RowIter, ColIter), FWHM_JumpMatrix(RowIter, ColIter) ] = ComputeAmpFWHM( normPDF, OffSet_thres, rads );
            % count the number of total angles: 
            NumAngles_JumpMatrix(RowIter, ColIter) = round(0.5*numel(MergedAngleMatrix{RowIter, ColIter}));           
        else 
            % fill in NaNs if there is not enough data
            Amp_JumpMatrix(RowIter, ColIter) = 0;
            FWHM_JumpMatrix(RowIter, ColIter) = 0;
            AC_JumpMatrix(RowIter, ColIter) = 0;
            f_180_0_JumpMatrix(RowIter, ColIter) = 0;
        end
    end
end












end

