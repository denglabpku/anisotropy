function [ output_struct ] = angleFWHM_Amp_HMM_analyzer_v5(input_struct )
%ANGLEFWHM_AMP_HMM_ANALYZER returns angle Amplitude and FWHM
%   written by Anders Sejr Hansen (AndersSejrHansen@post.harvard.edu;
%   @Anders_S_Hansen; https://anderssejrhansen.wordpress.com)
%   License: GNU GPL v3
%   This function takes the data and returns the angle Amplitude and FWHM
%   as well does subsampling on the data

%   UPDATE - Version 2; July 10, 2017
%   Enable analysis at multiple length scales

%   UPDATE - Version 3; July 27, 2017
%   Analyze each trajectory individually and the use QC-threshold to get
%   rid of weird outlier trajectories

%   UPDATE - Version 4: October 14, 2017
%   Updated for analysis of spatiotemporal analysis and also included the
%   fold-change (180 +/- 30 degrees) / (0 +/- 30 degrees) as output. This
%   is nearly identical to the Asymmetry Coefficient, but maybe a bit more
%   intuitive. 
%   Also added an upper limit on jumps that will be considered.
%   Added this upper limit to all angle calculations

%   UPDATE - Version 5: October 26, 2017
%   Updated to include the following:
%   -   count number of angles:
%       - for all jumps
%       - per displacement length bin
%   This will be helpful for taking weighted averages
%   - set global min-min threshold: only if jumps exceed this threshold
%   will they be counted
%   Also, record the full matrix of angles: 
%   Columns: 1st jump (length)
%   Rows:    2nd jump (length)
%   Colour:  Either Amp, f_180_0, FWHM


% define a structure for returning the outputs
output_struct = struct;

% conveniently store the inputs
CellTracks = input_struct(1).CellTracks;
CellTrackViterbiClass = input_struct(1).CellTrackViterbiClass;
BoundState = input_struct(1).BoundState;
FreeState = input_struct(1).FreeState;
nBins = input_struct(1).nBins;
rads = input_struct(1).rads;
MinMinJumpThres = input_struct(1).MinMinJumpThres;
GlobalMinJumpThres = input_struct(1).GlobalMinJumpThres;
minFree = input_struct(1).minFree;
OffSet_thres = input_struct(1).OffSet_thres;
JackKnife_fraction = input_struct(1).JackKnife_fraction;
JackKnife_iterations = input_struct(1).JackKnife_iterations;
minAngleNumber = input_struct(1).minAngleNumber;
MinNumAngles = input_struct(1).MinNumAngles;
MaxAsymAnglesFrac = input_struct(1).MaxAsymAnglesFrac;
MaxJump = input_struct(1).MaxJump;
MovingThreshold = input_struct(1).MovingThreshold;

% Asym conditions
minAsym = pi-1.5/18*pi; %165 degree
maxAsym = pi+1.5/18*pi; %195 degree
% create matrix to store all angles as a function of distances:

% MatAngles is an Nx4 matrix:
% column 1: first displacement length
% column 2: second displacement length
% column 3: first angle
% column 4: 2*pi-angle
MatAngles = [];



%% %%%%%%%%%%%%%%%% SLICE UP TRAJECTORIES AND COMPUTE ANGLES %%%%%%%%%%%%%%%%%
AnglesMinMinThres = [];
AllAngles = [];
%Find all trajectory pieces that have a length of at least 3 frames in the
%moving configuration. 
for i=1:size(CellTracks,2)
    %Check to see if there are at least two jumps classified as free in the
    %trajectory
    if length(find(CellTrackViterbiClass{i} == FreeState)) >= minFree
        %Now you need to make sure that the jumps are adjacent. 
        FreeIdx = find(CellTrackViterbiClass{i} == FreeState);
        ranges = cell(0); %make a cell array to store the vectors

        %Now comes the complicated part: there can be multiple switches
        %between the free and bound states. So you need to figure out the
        %indices of these switches to ensure that you are only counting
        %adjacent Free parts.
        
        %Zuhui Comments: for the "base" dataset that are not subsampled,
        %tracks containing bound and free mixture will also be considered
        %when calculate angle, but only limited to the free segments of the
        %track.
        
        %Figure out how many switches there are by how many non-adjacent
        %elements you have en Free Idx:
        StateSwitches = find(diff(FreeIdx)>1); %Number of switches
        if isempty(StateSwitches)
            %There were no switches, so the range is simply FreeIdx:
            ranges{1} = FreeIdx;
        else
            %There were switches: Now find out if there was only one or if
            %there were more
            if length(StateSwitches) == 1
                ranges{1} = FreeIdx(1):FreeIdx(StateSwitches);
                ranges{2} = FreeIdx(StateSwitches+1):FreeIdx(end);
            elseif length(StateSwitches) > 1
                %There are now at least three segments of free particles
                ranges{1} = FreeIdx(1:StateSwitches(1));
                for k=1:length(StateSwitches)-1
                    ranges{k+1} = FreeIdx(StateSwitches(k)+1):FreeIdx(StateSwitches(k+1));
                end
                ranges{length(ranges)+1} = FreeIdx(StateSwitches(end)+1):FreeIdx(end);
            end
        end

        temp_DistAngles = [];
        for k=1:length(ranges)
            %Now loop through all of the ranges of free particles and save
            %the angles
            CurrFreeIdx = ranges{k};

            %Now perform the angle calculations:
            for j=1:length(CurrFreeIdx)-1
                %define the x,y-coordinates of the three points
                p1 = [CellTracks{i}(CurrFreeIdx(j),1), CellTracks{i}(CurrFreeIdx(j),2)];
                p2 = [CellTracks{i}(CurrFreeIdx(j)+1,1), CellTracks{i}(CurrFreeIdx(j)+1,2)];
                p3 = [CellTracks{i}(CurrFreeIdx(j)+2,1), CellTracks{i}(CurrFreeIdx(j)+2,2)];

                %Check if the points meet the moving threshold
                Distances = [pdist([p1; p2]) pdist([p2; p3])];

                %now define the vectors as column vectors
                v1 = (p2-p1)';
                v2 = (p3-p2)';
                %Calculate the angle: make it symmetric
                angle = zeros(1,2);
                angle(1,1) = abs(atan2(det([v1,v2]),dot(v1,v2)));
                %this angle is in radians. Make it symmetric by adding
                %it
                angle(1,2) =  2*pi-angle(1,1);
                
                % save the stuff: 
                temp_DistAngles = vertcat(temp_DistAngles, [Distances(1), Distances(2), angle]);
                
                
                
            end
        end
        if ~isempty(temp_DistAngles)
            curr_angles = temp_DistAngles(:,3:4);
            curr_angles = curr_angles(:)';
            % calculate if the trajectory is weird
            NumbAsym = length(find(curr_angles>minAsym & curr_angles < maxAsym));
            FracAsym =  NumbAsym / length(curr_angles);
            % should we use the trajectory
            use_traj = 1;
            if NumbAsym > MinNumAngles
                if FracAsym > MaxAsymAnglesFrac
                    use_traj = 0;
                    %disp(['Trajectory number ', num2str(i), ' was excluded from analysis']);
                end
            end

            % test QC: only use if passes QC
            if use_traj == 1
                
                % test to see if each angle was bigger than the minimum
                % threshold
                for MinIter = 1:size(temp_DistAngles,1)
                    
                    %%%%%% THRESHOLD FOR AGGREGATE ANALYSIS OF ANGLES %%%%%
                    %ZW comments: below treshold also used in multitemporal
                    %scale analysis
                    %If the minimum jump is bigger than the Minimial
                    %MinThreshold then save the gave
                    if min(temp_DistAngles(MinIter,1:2)) > MinMinJumpThres
                        if max(temp_DistAngles(MinIter,1:2)) < MaxJump
                            AnglesMinMinThres = horzcat(AnglesMinMinThres, temp_DistAngles(MinIter,3:4));  
                        end
                    end
                    
                    %%%%%% GLOBAL THRESHOLD FOR SPATIOTEMPORAL ANALYSIS OF ANGLES %%%%%
                    %ZW comments: below treshold only used in multispatial
                    %scale analysis NOT in multitemporal scale analysis!!! 
                    % if jump length passed global minimal threshold:
                    if min(temp_DistAngles(MinIter,1:2)) > GlobalMinJumpThres
                        if max(temp_DistAngles(MinIter,1:2)) < MaxJump
                            
                            % store both jumps and angles:
                            MatAngles = vertcat(MatAngles, temp_DistAngles(MinIter,:));
                            %Store the angles but not jumps
                            AllAngles = horzcat(AllAngles, temp_DistAngles(MinIter,3:4));
                        end
                    end
                    
                    
                end

            end
        end
    end
end

% define key vectors:
jack_Amp = zeros(1,JackKnife_iterations);
jack_FWHM = zeros(1,JackKnife_iterations);
jack_AC = zeros(1,JackKnife_iterations);
jack_f_180_0 = zeros(1,JackKnife_iterations);

% what if there are no jumps, see else
%  ensure there are enough angle data, below this number, don't do calculations
if length(AnglesMinMinThres) > MinNumAngles

    %%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTE THE PDF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate the Assymmetry Coefficient for all of the data
    [Alla,AllVals] = rose(AnglesMinMinThres, nBins);
    AllVals = AllVals./sum(AllVals);
    %Calculate the Assymmetry Coefficient (Izeddin et al., eLife, 2014)
    %Use 0 +/- 30 degrees and 180 +/- 30 degrees as the 
    NonRedundantVals = AllVals(2:4:end);
    NumElements = round(length(NonRedundantVals)*30/360);
    For = [1:1:NumElements length(NonRedundantVals)-NumElements+1:1:length(NonRedundantVals)]; % add "+1" to get correct FOR array 
    Back = [round(length(NonRedundantVals)/2)-NumElements+1:1:round(length(NonRedundantVals)/2)+NumElements];
    SingleAsymCoef = log2(mean(NonRedundantVals(For))/mean(NonRedundantVals(Back)));
    Single_f_180_0 = mean(NonRedundantVals(Back)) / mean(NonRedundantVals(For));

    % make a normalised PDF:
    mean_normPDF = NonRedundantVals ./sum(NonRedundantVals);

    % save Alla:
    output_struct(1).Alla = Alla;
    output_struct(1).AllVals = AllVals;


    %%%%%%%%%%%%%%%%%%%%%%% COMPUTE AMPLITUDE AND FWHM %%%%%%%%%%%%%%%%%%%%%%%%
    [ mean_Amp, mean_FWHM_degree ] = ComputeAmpFWHM( mean_normPDF, OffSet_thres, rads );

    %%% JACK-KNIFE RE-SAMPLING TO ESTIMATE ERROR %%%
    FracToSample = round( length(AnglesMinMinThres) * JackKnife_fraction );
    for iter = 1:JackKnife_iterations
        % take a sub-sample
        subsampled_minThresAngles = datasample(AnglesMinMinThres, FracToSample, 'Replace', false);

        clear Alla AllVals NonRedundantVals
        % calculate normPDF:
        [Alla,AllVals] = rose(subsampled_minThresAngles, nBins);
        AllVals = AllVals./sum(AllVals);
        NonRedundantVals = AllVals(2:4:end);
        % make a normalised PDF:
        subsampled_normPDF = NonRedundantVals ./sum(NonRedundantVals);

        % calculate Amplitude and FWHM
        [jack_Amp(1,iter), jack_FWHM(1,iter)] = ComputeAmpFWHM( subsampled_normPDF, OffSet_thres, rads );

        % calculate sub-sampled AC
        NumElements = round(length(NonRedundantVals)*30/360);
        For = [1:1:NumElements length(NonRedundantVals)-NumElements+1:1:length(NonRedundantVals)]; % % add "+1" to get correct FOR array
        Back = [round(length(NonRedundantVals)/2)-NumElements+1:1:round(length(NonRedundantVals)/2)+NumElements];
        jack_AC(1,iter) = log2(mean(NonRedundantVals(For))/mean(NonRedundantVals(Back)));
        jack_f_180_0(1,iter) =  mean(NonRedundantVals(Back)) / mean(NonRedundantVals(For));
    end
else
    % if there are no jumps, just hand back empty stuff:
    Alla = [];
    AllVals = [];
    output_struct(1).Alla = Alla;
    output_struct(1).AllVals = AllVals;
    mean_normPDF = [];
    mean_Amp = 0;
    mean_FWHM_degree = 0;
    SingleAsymCoef = 0;
    Single_f_180_0 = 0;
end
%%%%%%%%%%%%%%%%%% ALL DONE: COMPILE THE RELEVANT OUTPUT %%%%%%%%%%%%%%%%%%
output_struct(1).normPDF = mean_normPDF;
output_struct(1).mean_Amp = mean_Amp;
output_struct(1).mean_FWHM = mean_FWHM_degree;
output_struct(1).jack_Amp = jack_Amp;
output_struct(1).jack_FWHM = jack_FWHM;
output_struct(1).AC = SingleAsymCoef;
output_struct(1).jack_AC = jack_AC;
output_struct(1).f_180_0 = Single_f_180_0;
output_struct(1).jack_f_180_0 = jack_f_180_0;

% save the total number of angles:
output_struct(1).TotalNumAngles = size(MatAngles,1);
output_struct(1).TotalMinMinAngles = round(0.5*length(AnglesMinMinThres));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% PERFORM ANALYSIS AS A FUNCTION OF LENGTH SCALE %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% perform two types of analysis:
%   - based on the mean translocation
%   - based on the minimum translocation
% If MovingThreshold is 0.2:1; so first element should be 200-400 nm and last
% element should be more than 1 um
AngleAtClosestMeanJump = cell(1,length(MovingThreshold));
AngleAtClosestMinJump = cell(1,length(MovingThreshold));
AngleMatrix = cell(length(MovingThreshold),length(MovingThreshold));
% add abitrary large value to MovingThreshol for indexing:
tempMovingThreshold = [MovingThreshold, 10000];

% now loop through all angles:
for iter = 1:size(MatAngles,1)
    curr_mean = mean(MatAngles(iter,1:2));
    curr_min = min(MatAngles(iter,1:2));
    curr_max = max(MatAngles(iter,1:2));
    curr_angles = MatAngles(iter,3:4);
    
       
    % test to see which they are bigger than:
    if curr_max < MaxJump
        ColumnIdx = 0;
        RowIdx = 0;
        for n=2:length(tempMovingThreshold)
            % calculate vectors
            if curr_mean > tempMovingThreshold(1,n-1) && curr_mean < tempMovingThreshold(1,n)
                AngleAtClosestMeanJump{n-1} = horzcat(AngleAtClosestMeanJump{n-1}, curr_angles);
            end
            if curr_min > tempMovingThreshold(1,n-1) && curr_min < tempMovingThreshold(1,n)
                AngleAtClosestMinJump{n-1} = horzcat(AngleAtClosestMinJump{n-1}, curr_angles);
            end
            % calculate full matrix:
            if MatAngles(iter,1) > tempMovingThreshold(1,n-1) && MatAngles(iter,1) < tempMovingThreshold(1,n)
                ColumnIdx = n-1;
            end
            if MatAngles(iter,2) > tempMovingThreshold(1,n-1) && MatAngles(iter,2) < tempMovingThreshold(1,n)
                RowIdx = n-1;
            end
            
        end
        % now see if you found both a column and row Idx:
        if ColumnIdx > 0 && RowIdx > 0
            AngleMatrix{ColumnIdx, RowIdx} = horzcat(AngleMatrix{ColumnIdx, RowIdx}, curr_angles);
        end
    end
end

% calculate Amp, AC and FWHM as a function of length
AmpAtClosestMean = zeros(1, length(MovingThreshold));
FWHMAtClosestMean = zeros(1, length(MovingThreshold));
ACAtClosestMean = zeros(1, length(MovingThreshold));
f_180_0_AtClosestMean = zeros(1, length(MovingThreshold));
AmpAtClosestMin = zeros(1, length(MovingThreshold));
FWHMAtClosestMin = zeros(1, length(MovingThreshold));
ACAtClosestMin = zeros(1, length(MovingThreshold));
f_180_0_AtClosestMin = zeros(1, length(MovingThreshold));
NumAnglesAtClosestMean = zeros(1, length(MovingThreshold));
NumAnglesAtClosestMin = zeros(1, length(MovingThreshold));
% for Jack-Knife
jack_AmpAtClosestMean = zeros(JackKnife_iterations, length(MovingThreshold));
jack_FWHMAtClosestMean = zeros(JackKnife_iterations, length(MovingThreshold));
jack_ACAtClosestMean = zeros(JackKnife_iterations, length(MovingThreshold));
jack_f_180_0_AtClosestMean = zeros(JackKnife_iterations, length(MovingThreshold));
jack_AmpAtClosestMin = zeros(JackKnife_iterations, length(MovingThreshold));
jack_FWHMAtClosestMin = zeros(JackKnife_iterations, length(MovingThreshold));
jack_ACAtClosestMin = zeros(JackKnife_iterations, length(MovingThreshold));
jack_f_180_0_AtClosestMin = zeros(JackKnife_iterations, length(MovingThreshold));

% Calculate Amp, AC and FWHM at Closest Mean and Min:
for iter = 1:length(MovingThreshold)
 
    %%%%%%%%%%%% CALCULATION FOR CLOSEST MEAN %%%%%%%%%%%%
    % ensure that there is enough data
    if numel(AngleAtClosestMeanJump{iter}) > minAngleNumber
        [Alla,AllVals] = rose(AngleAtClosestMeanJump{iter}, nBins);
        AllVals = AllVals./sum(AllVals);
        %Calculate the Assymmetry Coefficient (Izeddin et al.)
        %Use 0 +/- 30 degrees and 180 +/- 30 degrees as the 
        NonRedundantVals = AllVals(2:4:end);
        NumElements = round(length(NonRedundantVals)*30/360);
        For = [1:1:NumElements length(NonRedundantVals)-NumElements+1:1:length(NonRedundantVals)]; % add "+1" to get correct FOR array
        Back = [round(length(NonRedundantVals)/2)-NumElements+1:1:round(length(NonRedundantVals)/2)+NumElements];
        ACAtClosestMean(1,iter) = log2(mean(NonRedundantVals(For))/mean(NonRedundantVals(Back)));
        f_180_0_AtClosestMean(1,iter) = mean(NonRedundantVals(Back))/mean(NonRedundantVals(For));
        % make a normalised PDF:
        normPDF = NonRedundantVals ./sum(NonRedundantVals);
        [ AmpAtClosestMean(1,iter), FWHMAtClosestMean(1,iter) ] = ComputeAmpFWHM( normPDF, OffSet_thres, rads );
        % count the number of total angles: 
        NumAnglesAtClosestMean(1,iter) = numel(AngleAtClosestMeanJump{iter});
    else
        % fill in NaNs if there is not enough data
        ACAtClosestMean(1,iter) = NaN;
        f_180_0_AtClosestMean(1,iter) = NaN;
        AmpAtClosestMean(1,iter) = NaN;
        FWHMAtClosestMean(1,iter) = NaN;
    end

    
    %%%%%%%%%%%% CALCULATION FOR CLOSEST MIN %%%%%%%%%%%%
    if numel(AngleAtClosestMinJump{iter}) > minAngleNumber
        [Alla,AllVals] = rose(AngleAtClosestMinJump{iter}, nBins);
        AllVals = AllVals./sum(AllVals);
        %Calculate the Assymmetry Coefficient (Izeddin et al.)
        %Use 0 +/- 30 degrees and 180 +/- 30 degrees as the 
        NonRedundantVals = AllVals(2:4:end);
        NumElements = round(length(NonRedundantVals)*30/360);
        For = [1:1:NumElements length(NonRedundantVals)-NumElements+1:1:length(NonRedundantVals)]; % add "+1" to get correct FOR array
        Back = [round(length(NonRedundantVals)/2)-NumElements+1:1:round(length(NonRedundantVals)/2)+NumElements];
        ACAtClosestMin(1,iter) = log2(mean(NonRedundantVals(For))/mean(NonRedundantVals(Back)));
        f_180_0_AtClosestMin(1,iter) = mean(NonRedundantVals(Back))/mean(NonRedundantVals(For));
        % make a normalised PDF:
        normPDF = NonRedundantVals ./sum(NonRedundantVals);
        [ AmpAtClosestMin(1,iter), FWHMAtClosestMin(1,iter) ] = ComputeAmpFWHM( normPDF, OffSet_thres, rads );
        % count the number of total angles: 
        NumAnglesAtClosestMin(1,iter) = numel(AngleAtClosestMinJump{iter});
    else
        % fill in NaNs if there is not enough data
        ACAtClosestMin(1,iter) = NaN;
        f_180_0_AtClosestMin(1,iter) = NaN;
        AmpAtClosestMin(1,iter) = NaN;
        FWHMAtClosestMin(1,iter) = NaN;
    end
    
    
    
    %%%%%%% NOW JACK-KNIFE SUB-SAMPLE FOR THE CURRENT LENGTH SCALE %%%%%%%
    FracToSample_mean = round( length(AngleAtClosestMeanJump{iter}) * JackKnife_fraction );
    FracToSample_min = round( length(AngleAtClosestMinJump{iter}) * JackKnife_fraction );
    for JackIter = 1:JackKnife_iterations
        
        %%%%%%%%%%%% CALCULATION FOR CLOSEST MEAN %%%%%%%%%%%%
        % take a sub-sample
        subsampled_minThresAngles = datasample(AngleAtClosestMeanJump{iter}, FracToSample_mean, 'Replace', false);
        % make sure that there is enough data
        if numel(subsampled_minThresAngles) > minAngleNumber
            % calculate normPDF:
            [Alla,AllVals] = rose(subsampled_minThresAngles, nBins);
            AllVals = AllVals./sum(AllVals);
            NonRedundantVals = AllVals(2:4:end);
            % make a normalised PDF:
            subsampled_normPDF = NonRedundantVals ./sum(NonRedundantVals);

            % calculate Amplitude and FWHM
            [jack_AmpAtClosestMean(JackIter,iter), jack_FWHMAtClosestMean(JackIter,iter)] = ComputeAmpFWHM( subsampled_normPDF, OffSet_thres, rads );

            % calculate sub-sampled AC
            NumElements = round(length(NonRedundantVals)*30/360);
            For = [1:1:NumElements length(NonRedundantVals)-NumElements+1:1:length(NonRedundantVals)]; % add "+1" to get correct FOR array
            Back = [round(length(NonRedundantVals)/2)-NumElements+1:1:round(length(NonRedundantVals)/2)+NumElements];
            jack_ACAtClosestMean(JackIter,iter) = log2(mean(NonRedundantVals(For))/mean(NonRedundantVals(Back)));
            jack_f_180_0_AtClosestMean(JackIter,iter) = mean(NonRedundantVals(Back))/mean(NonRedundantVals(For));
        else
            % fill in NaNs if there is not enough data
            jack_ACAtClosestMean(JackIter,iter) = NaN;
            jack_f_180_0_AtClosestMean(JackIter,iter) = NaN;
            jack_AmpAtClosestMean(JackIter,iter) = NaN;
            jack_FWHMAtClosestMean(JackIter,iter) = NaN;
        end
        
        %%%%%%%%%%%% CALCULATION FOR CLOSEST MIN %%%%%%%%%%%%
        % take a sub-sample
        subsampled_minThresAngles = datasample(AngleAtClosestMinJump{iter}, FracToSample_min, 'Replace', false);
        % make sure that there is enough data
        if numel(subsampled_minThresAngles) > minAngleNumber
            % calculate normPDF:
            [Alla,AllVals] = rose(subsampled_minThresAngles, nBins);
            AllVals = AllVals./sum(AllVals);
            NonRedundantVals = AllVals(2:4:end);
            % make a normalised PDF:
            subsampled_normPDF = NonRedundantVals ./sum(NonRedundantVals);

            % calculate Amplitude and FWHM
            [jack_AmpAtClosestMin(JackIter,iter), jack_FWHMAtClosestMin(JackIter,iter)] = ComputeAmpFWHM( subsampled_normPDF, OffSet_thres, rads );

            % calculate sub-sampled AC
            NumElements = round(length(NonRedundantVals)*30/360);
            For = [1:1:NumElements length(NonRedundantVals)-NumElements+1:1:length(NonRedundantVals)]; % add "+1" to get correct FOR array
            Back = [round(length(NonRedundantVals)/2)-NumElements+1:1:round(length(NonRedundantVals)/2)+NumElements];
            jack_ACAtClosestMin(JackIter,iter) = log2(mean(NonRedundantVals(For))/mean(NonRedundantVals(Back)));
            jack_f_180_0_AtClosestMin(JackIter,iter) = mean(NonRedundantVals(Back))/mean(NonRedundantVals(For));
        else
            % fill in NaNs if there is not enough data
            jack_f_180_0_AtClosestMin(JackIter,iter) = NaN;
            jack_ACAtClosestMin(JackIter,iter) = NaN;
            jack_AmpAtClosestMin(JackIter,iter) = NaN;
            jack_FWHMAtClosestMin(JackIter,iter) = NaN;
        end
    end
end




%%%%%%%%%%%%%%%%%% ALL DONE: COMPILE THE RELEVANT OUTPUT %%%%%%%%%%%%%%%%%%

% at multiple length scales
output_struct(1).AmpAtClosestMean = AmpAtClosestMean;
output_struct(1).FWHMAtClosestMean = FWHMAtClosestMean;
output_struct(1).ACAtClosestMean = ACAtClosestMean;
output_struct(1).f_180_0_AtClosestMean = f_180_0_AtClosestMean;
output_struct(1).AmpAtClosestMin = AmpAtClosestMin;
output_struct(1).FWHMAtClosestMin = FWHMAtClosestMin;
output_struct(1).ACAtClosestMin = ACAtClosestMin;
output_struct(1).f_180_0_AtClosestMin = f_180_0_AtClosestMin;
output_struct(1).AngleMatrix = AngleMatrix;

% total number of angles for doing weighted calculations:
output_struct(1).NumAnglesAtClosestMean = NumAnglesAtClosestMean;
output_struct(1).NumAnglesAtClosestMin = NumAnglesAtClosestMin;

% jack-knife sampled
output_struct(1).jack_AmpAtClosestMean = jack_AmpAtClosestMean;
output_struct(1).jack_FWHMAtClosestMean = jack_FWHMAtClosestMean;
output_struct(1).jack_ACAtClosestMean = jack_ACAtClosestMean;
output_struct(1).jack_f_180_0_AtClosestMean = jack_f_180_0_AtClosestMean;
output_struct(1).jack_AmpAtClosestMin = jack_AmpAtClosestMin;
output_struct(1).jack_FWHMAtClosestMin = jack_FWHMAtClosestMin;
output_struct(1).jack_ACAtClosestMin = jack_ACAtClosestMin;
output_struct(1).jack_f_180_0_AtClosestMin = jack_f_180_0_AtClosestMin;


end

