function outputStruct = FreeSegmentAsymCalculation(inputStruct)
%FreeSegmentAsymCalculation Find the defined number (FreSeg_size) of free segment before and after a bound segment and calculate the corresponding f_180_0.
%   Detailed explanation goes here:
%   (1)Determine which track are mixed with bound and unbound segments
%   (2)Find the free segment after or before a bound state segment
%   (3)Calculate angle of consecutive free segments 
%   (4)Calculate f_180_0 of free segments after bound segment
%   (5)Calculate f_180_0 of free segments before bound segment
%   (6)Calculate f_180_0 of pooled free segments after and before bound segment


%   Zuhui Wang 2021/09/17

%% INTENTIONALLY LEFT BLANK
%%
load([inputStruct.dataPath inputStruct.dataName],'CellTracks','CellTrackViterbiClass');
fprintf('Data loaded:\n %s%s\n.',inputStruct.dataPath, inputStruct.dataName);

BoundState = inputStruct.BoundState;
FreeState = inputStruct.FreeState;
FreSeg_size = inputStruct.FreSeg_size;
nBins = inputStruct.nBins;
JackKnife_fraction = inputStruct.JackKnife_fraction;
JackKnife_iterations = inputStruct.JackKnife_iterations;
MinMinJumpThres = inputStruct.MinMinJumpThres;
MaxJump = inputStruct.MaxJump;

%% Determine which track are mixed with bound and unbound segments
BoUnboMixTrackIdx = cellfun(@CheckUnique,CellTrackViterbiClass);
BoUnboMixTrackIdx = find(BoUnboMixTrackIdx);

BoUnbo_CellTrackViterbiClass = CellTrackViterbiClass(BoUnboMixTrackIdx);

BoUnbo_CellTracks = CellTracks(BoUnboMixTrackIdx);


%% Find the free segment after or before a bound state segment

j = 1; % iter of index of free segment in freesegment after bound segment
k = 1; % iter of index of free segment in freesegment before bound segment

FreeSeg_AfBo = {}; % the index of Free segment after bound state segment, same as 'BoUnbo_CellTrackViterbiClass'
FreeSeg_BfBo = {}; % the index of Free segment before bound state segment
BoSegAfterFreeSeg_BfB = {}; % Consecutive bound segments after the free segments (those defined in FreeSeg_BfBo)
IdxLUT_FreeSeg_AfBo = []; % the same length as 'FreeSeg_AfBo', but tell you the index of track index in 'BoUnbo_CellTrackViterbiClass'
IdxLUT_FreeSeg_BfBo = []; % the same length as 'FreeSeg_BfBo', but tell you the index of track index in 'BoUnbo_CellTrackViterbiClass'

totBoUnboNum = length(BoUnbo_CellTracks);

for i = 1:totBoUnboNum
    
    gap = diff(double(BoUnbo_CellTrackViterbiClass{i}));

    % Find the index of consecutive free segment after a bound state
    if sum(gap == 1) % if Bound to free transition exist
        BoFreIdx = find(gap == 1);
        for BoFreIdx_iter = BoFreIdx'
            seg1_fre_idx = BoFreIdx_iter + 1;
            seg2_fre_idx = BoFreIdx_iter + 2;
            FreeSeg_AfBo{j} = [seg1_fre_idx];
            if seg2_fre_idx > length(BoUnbo_CellTrackViterbiClass{i})
                continue;
            else
                while BoUnbo_CellTrackViterbiClass{i}(seg2_fre_idx) == FreeState
                    if length(FreeSeg_AfBo{j}) >= FreSeg_size % exceed defined max free segment size
                        break;
                    else
                        FreeSeg_AfBo{j} = [FreeSeg_AfBo{j}  seg2_fre_idx];
                        seg2_fre_idx = seg2_fre_idx + 1; % next segment
                        if seg2_fre_idx > length(BoUnbo_CellTrackViterbiClass{i})%NOT exceed the max length of state array
                           break;
                        end
                    end
                end 
                IdxLUT_FreeSeg_AfBo = [IdxLUT_FreeSeg_AfBo,i]; % register the current Track idx of FreeSeg_AfBo
                j = j + 1;
            end
        end
    end

    
    % Find the index of consecutive free segment before a bound state
    if sum(gap == -1)
        BoFreIdx = find(gap == -1);
        for BoFreIdx_iter = BoFreIdx'
            seg1_fre_idx = BoFreIdx_iter;
            segM1_fre_idx = BoFreIdx_iter - 1; % segM1_fre_idx means seg minus 1 free segment, i.e. segment before the seg1_fre, which is the last free segment before bound segment
            FreeSeg_BfBo{k} = [seg1_fre_idx];
            
            segP1_fre_idx = BoFreIdx_iter + 1;
            BoSegAfterFreeSeg_BfB{k} = [segP1_fre_idx];
              
            if segM1_fre_idx < 1 % NOT exceed first element of state array
                continue;
                %   We do not include the situation when only one free segment exists before
                %   the bound segment. 
            else
                % Continue add free segments if the segment before the last
                % added free segment are still free and not exceed the
                % first segments of the state array
                while BoUnbo_CellTrackViterbiClass{i}(segM1_fre_idx) == FreeState
                    if length(FreeSeg_BfBo{k}) >= FreSeg_size % exceed defined max free segment size
                        break;
                    else 
                        FreeSeg_BfBo{k} = [segM1_fre_idx FreeSeg_BfBo{k}];
                        segM1_fre_idx = segM1_fre_idx - 1;
                        if segM1_fre_idx < 1 % NOT exceed first element of state array
                            break;
                        end
                    end
                end
                
                while (segP1_fre_idx+1)<= length(BoUnbo_CellTrackViterbiClass{i})
                    if BoUnbo_CellTrackViterbiClass{i}(segP1_fre_idx+1) == BoundState
                        BoSegAfterFreeSeg_BfB{k} = [BoSegAfterFreeSeg_BfB{k} segP1_fre_idx+1];
                        segP1_fre_idx = segP1_fre_idx +1;
                    else
                        break
                    end
                    
                end
                
                IdxLUT_FreeSeg_BfBo = [IdxLUT_FreeSeg_BfBo,i]; % register the current Track idx of FreeSeg_BfBo
                k = k + 1;
            end
        end
    end
    
    fprintf('The free segments after and before a bound segment have been split from bound-unbound mixed tracks: %d / %d.\n',i,totBoUnboNum);
 
end

% Remove the recording that contain only ONE free segment
temp_free_length_AfBo = cellfun(@length,FreeSeg_AfBo);
filt_FreeSeg_AfBo = FreeSeg_AfBo(temp_free_length_AfBo ~= 1);
filt_IdxLUT_FreeSeg_AfBo = IdxLUT_FreeSeg_AfBo(temp_free_length_AfBo ~= 1);

temp_free_length_BfBo = cellfun(@length,FreeSeg_BfBo);
filt_FreeSeg_BfBo = FreeSeg_BfBo(temp_free_length_BfBo ~= 1);
filt_IdxLUT_FreeSeg_BfBo = IdxLUT_FreeSeg_BfBo(temp_free_length_BfBo ~= 1);
filt_BoSegAfterFreeSeg_BfB = BoSegAfterFreeSeg_BfB(temp_free_length_BfBo ~= 1);

%% Calculate angle of consecutive segments
%%%============================================%%%%%%%%%%
% Calculate the angle of free segments after bound segment
k = 1;
totFreSeg_AfBo = length(filt_FreeSeg_AfBo);
for i = 1:length(filt_IdxLUT_FreeSeg_AfBo)
    if length(filt_FreeSeg_AfBo{i}) == FreSeg_size % only calculate angles from FreSeg_size number of free segments
     
        trackStartIdx = filt_FreeSeg_AfBo{i}(1);
        trackEndIdx = filt_FreeSeg_AfBo{i}(end);

        for j = 1:length(filt_FreeSeg_AfBo{i})-1

            p1 = [BoUnbo_CellTracks{filt_IdxLUT_FreeSeg_AfBo(i)}(trackStartIdx+j-1,1),BoUnbo_CellTracks{filt_IdxLUT_FreeSeg_AfBo(i)}(trackStartIdx+j-1,2)];
            p2 = [BoUnbo_CellTracks{filt_IdxLUT_FreeSeg_AfBo(i)}(trackStartIdx+j-1+1,1),BoUnbo_CellTracks{filt_IdxLUT_FreeSeg_AfBo(i)}(trackStartIdx+j-1+1,2)];
            p3 = [BoUnbo_CellTracks{filt_IdxLUT_FreeSeg_AfBo(i)}(trackStartIdx+j-1+2,1),BoUnbo_CellTracks{filt_IdxLUT_FreeSeg_AfBo(i)}(trackStartIdx+j-1+2,2)];

            %Calculate two consecutive jump lengths
            Distances = nan(1,2);
            Distances(1:2) = [pdist([p1; p2]) pdist([p2; p3])];  

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
            % [ BoUnboTrack index, first jump length, second jump
            % length, angle in Rad, complement angle in Rad]
            if min(Distances) > MinMinJumpThres
                if max(Distances) < MaxJump
                   angleMatrix_FreAfBo(k,:) = [filt_IdxLUT_FreeSeg_AfBo(i), Distances, angle];
                   ffilt_FreeSeg_AfBo(k) = filt_FreeSeg_AfBo(i);
                   ffilt_IdxLUT_FreeSeg_AfBo(k) = filt_IdxLUT_FreeSeg_AfBo(i);
                   k = k+1;
                end
            end
        end
    end
    fprintf('Current processing %d / %d spilt free segments after bound state segments.\n',i,totFreSeg_AfBo);
end

%%%============================================%%%%%%%%%%
% Calculate the angle of free segments before bound segment
k = 1;
totFreSeg_BfBo = length(filt_FreeSeg_BfBo);
for i = 1:length(filt_IdxLUT_FreeSeg_BfBo)
    if length(filt_FreeSeg_BfBo{i}) == FreSeg_size % only calculate angles from FreSeg_size number of free segments
     
        trackStartIdx = filt_FreeSeg_BfBo{i}(1);
        trackEndIdx = filt_FreeSeg_BfBo{i}(end);

        for j = 1:length(filt_FreeSeg_BfBo{i})-1

            p1 = [BoUnbo_CellTracks{filt_IdxLUT_FreeSeg_BfBo(i)}(trackStartIdx+j-1,1),BoUnbo_CellTracks{filt_IdxLUT_FreeSeg_BfBo(i)}(trackStartIdx+j-1,2)];
            p2 = [BoUnbo_CellTracks{filt_IdxLUT_FreeSeg_BfBo(i)}(trackStartIdx+j-1+1,1),BoUnbo_CellTracks{filt_IdxLUT_FreeSeg_BfBo(i)}(trackStartIdx+j-1+1,2)];
            p3 = [BoUnbo_CellTracks{filt_IdxLUT_FreeSeg_BfBo(i)}(trackStartIdx+j-1+2,1),BoUnbo_CellTracks{filt_IdxLUT_FreeSeg_BfBo(i)}(trackStartIdx+j-1+2,2)];

            %Calculate two consecutive jump lengths
            Distances = nan(1,2);
            Distances(1:2) = [pdist([p1; p2]) pdist([p2; p3])];  

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
            % [ BoUnboTrack index, first jump length, second jump
            % length, angle in Rad, complement angle in Rad]
            if min(Distances) > MinMinJumpThres
                if max(Distances) < MaxJump
                   angleMatrix_FreBfBo(k,:) = [filt_IdxLUT_FreeSeg_BfBo(i), Distances, angle];
                   ffilt_FreeSeg_BfBo(k) = filt_FreeSeg_BfBo(i);
                   ffilt_IdxLUT_FreeSeg_BfBo(k) = filt_IdxLUT_FreeSeg_BfBo(i);
                   ffilt_BoSegAfterFreeSeg_BfB(k) = filt_BoSegAfterFreeSeg_BfB(i);
                   k = k+1;
                end
            end
        end
    end
    fprintf('Current processing %d / %d spilt free segments before bound state segments.\n',i,totFreSeg_BfBo);
end


%% Calculate f_180_0 of free segments after bound segment

% nBins = 1*36; %so 10 degrees per bin
AnglesMinMinThres = [angleMatrix_FreAfBo(:,4);angleMatrix_FreAfBo(:,5)];

%%%%%%%%%%%%%%%%%% Calculate f_180_0 %%%%%%%%%%%%%%%%%%%%%
%Calculate the Assymmetry Coefficient for all of the data
[~,AllVals] = rose(AnglesMinMinThres, nBins);
AllVals = AllVals./sum(AllVals);
%Calculate the Assymmetry Coefficient (Izeddin et al., eLife, 2014)
%Use 0 +/- 30 degrees and 180 +/- 30 degrees as the 
NonRedundantVals = AllVals(2:4:end);
NumElements = round(length(NonRedundantVals)*30/360);
For = [1:1:NumElements length(NonRedundantVals)-NumElements+1:1:length(NonRedundantVals)]; % add "+1" to get correct FOR array 
Back = [round(length(NonRedundantVals)/2)-NumElements+1:1:round(length(NonRedundantVals)/2)+NumElements];
% SingleAsymCoef = log2(mean(NonRedundantVals(For))/mean(NonRedundantVals(Back))); % the assym coefficient used in Izeddin et al, 2014.
FreAfBo_Single_f_180_0 = mean(NonRedundantVals(Back)) / mean(NonRedundantVals(For));


%%%%%%%%%%%%%%%%%% bootstrap to obtain std of f_180_0 %%%%%%%%%%%%%%%%%%%%%
% % JackKnife params
% JackKnife_fraction = 0.5;
% JackKnife_iterations = 50;
%%% JACK-KNIFE RE-SAMPLING TO ESTIMATE ERROR %%%
FracToSample = round( length(AnglesMinMinThres) * JackKnife_fraction );
for iter = 1:JackKnife_iterations
    % take a sub-sample
    subsampled_minThresAngles = datasample(AnglesMinMinThres, FracToSample, 'Replace', false);

    clear Alla AllVals NonRedundantVals
    % calculate normPDF:
    [~,AllVals] = rose(subsampled_minThresAngles, nBins);
    AllVals = AllVals./sum(AllVals);
    NonRedundantVals = AllVals(2:4:end);
    % make a normalised PDF:
%     subsampled_normPDF = NonRedundantVals ./sum(NonRedundantVals);

    % calculate sub-sampled AC
    NumElements = round(length(NonRedundantVals)*30/360);
    For = [1:1:NumElements length(NonRedundantVals)-NumElements+1:1:length(NonRedundantVals)]; % % add "+1" to get correct FOR array
    Back = [round(length(NonRedundantVals)/2)-NumElements+1:1:round(length(NonRedundantVals)/2)+NumElements];

    jack_f_180_0(1,iter) =  mean(NonRedundantVals(Back)) / mean(NonRedundantVals(For));
end

FreAfBo_std_f_180_0 = std(jack_f_180_0);

clear jack_f_180_0 For Back AnglesMinMinThres



%% Calculate f_180_0 of free segments before bound segment

% nBins = 1*36; %so 10 degrees per bin
AnglesMinMinThres = [angleMatrix_FreBfBo(:,4);angleMatrix_FreBfBo(:,5)];

%%%%%%%%%%%%%%%%%% Calculate f_180_0 %%%%%%%%%%%%%%%%%%%%%
%Calculate the Assymmetry Coefficient for all of the data
[~,AllVals] = rose(AnglesMinMinThres, nBins);
AllVals = AllVals./sum(AllVals);
%Calculate the Assymmetry Coefficient (Izeddin et al., eLife, 2014)
%Use 0 +/- 30 degrees and 180 +/- 30 degrees as the 
NonRedundantVals = AllVals(2:4:end);
NumElements = round(length(NonRedundantVals)*30/360);
For = [1:1:NumElements length(NonRedundantVals)-NumElements+1:1:length(NonRedundantVals)]; % add "+1" to get correct FOR array 
Back = [round(length(NonRedundantVals)/2)-NumElements+1:1:round(length(NonRedundantVals)/2)+NumElements];
% SingleAsymCoef = log2(mean(NonRedundantVals(For))/mean(NonRedundantVals(Back))); % the assym coefficient used in Izeddin et al, 2014.
FreBfBo_Single_f_180_0 = mean(NonRedundantVals(Back)) / mean(NonRedundantVals(For));


%%%%%%%%%%%%%%%%%% bootstrap to obtain std of f_180_0 %%%%%%%%%%%%%%%%%%%%%
% % JackKnife params
% JackKnife_fraction = 0.5;
% JackKnife_iterations = 50;
%%% JACK-KNIFE RE-SAMPLING TO ESTIMATE ERROR %%%
FracToSample = round( length(AnglesMinMinThres) * JackKnife_fraction );
for iter = 1:JackKnife_iterations
    % take a sub-sample
    subsampled_minThresAngles = datasample(AnglesMinMinThres, FracToSample, 'Replace', false);

    clear Alla AllVals NonRedundantVals
    % calculate normPDF:
    [~,AllVals] = rose(subsampled_minThresAngles, nBins);
    AllVals = AllVals./sum(AllVals);
    NonRedundantVals = AllVals(2:4:end);
    % make a normalised PDF:
%     subsampled_normPDF = NonRedundantVals ./sum(NonRedundantVals);

    % calculate sub-sampled AC
    NumElements = round(length(NonRedundantVals)*30/360);
    For = [1:1:NumElements length(NonRedundantVals)-NumElements+1:1:length(NonRedundantVals)]; % % add "+1" to get correct FOR array
    Back = [round(length(NonRedundantVals)/2)-NumElements+1:1:round(length(NonRedundantVals)/2)+NumElements];

    jack_f_180_0(1,iter) =  mean(NonRedundantVals(Back)) / mean(NonRedundantVals(For));
end

FreBfBo_std_f_180_0 = std(jack_f_180_0);

clear jack_f_180_0 For Back AnglesMinMinThres

%% Calculate f_180_0 of pooled free segments after and before bound segment
% nBins = 1*36; %so 10 degrees per bin
AnglesMinMinThres = [angleMatrix_FreAfBo(:,4);angleMatrix_FreAfBo(:,5);angleMatrix_FreBfBo(:,4);angleMatrix_FreBfBo(:,5)];

%%%%%%%%%%%%%%%%%% Calculate f_180_0 %%%%%%%%%%%%%%%%%%%%%
%Calculate the Assymmetry Coefficient for all of the data
[Alla,AllVals] = rose(AnglesMinMinThres, nBins);
AllVals = AllVals./sum(AllVals);
%Calculate the Assymmetry Coefficient (Izeddin et al., eLife, 2014)
%Use 0 +/- 30 degrees and 180 +/- 30 degrees as the 
NonRedundantVals = AllVals(2:4:end);
NumElements = round(length(NonRedundantVals)*30/360);
For = [1:1:NumElements length(NonRedundantVals)-NumElements+1:1:length(NonRedundantVals)]; % add "+1" to get correct FOR array 
Back = [round(length(NonRedundantVals)/2)-NumElements+1:1:round(length(NonRedundantVals)/2)+NumElements];
% SingleAsymCoef = log2(mean(NonRedundantVals(For))/mean(NonRedundantVals(Back))); % the assym coefficient used in Izeddin et al, 2014.
FreABfBo_Single_f_180_0 = mean(NonRedundantVals(Back)) / mean(NonRedundantVals(For));

% save Alla:
outputStruct(1).Alla = Alla; % theta(angle) for polarplot from all angles from pooled free segments after and before bound segment
outputStruct(1).AllVals = AllVals; % radius for polarplot from all angles from pooled free segments after and before bound segment


%%%%%%%%%%%%%%%%%% bootstrap to obtain std of f_180_0 %%%%%%%%%%%%%%%%%%%%%
% % JackKnife params
% JackKnife_fraction = 0.5;
% JackKnife_iterations = 50;
%%% JACK-KNIFE RE-SAMPLING TO ESTIMATE ERROR %%%
FracToSample = round( length(AnglesMinMinThres) * JackKnife_fraction );
for iter = 1:JackKnife_iterations
    % take a sub-sample
    subsampled_minThresAngles = datasample(AnglesMinMinThres, FracToSample, 'Replace', false);

    clear Alla AllVals NonRedundantVals
    % calculate normPDF:
    [~,AllVals] = rose(subsampled_minThresAngles, nBins);
    AllVals = AllVals./sum(AllVals);
    NonRedundantVals = AllVals(2:4:end);
    % make a normalised PDF:
%     subsampled_normPDF = NonRedundantVals ./sum(NonRedundantVals);

    % calculate sub-sampled AC
    NumElements = round(length(NonRedundantVals)*30/360);
    For = [1:1:NumElements length(NonRedundantVals)-NumElements+1:1:length(NonRedundantVals)]; % % add "+1" to get correct FOR array
    Back = [round(length(NonRedundantVals)/2)-NumElements+1:1:round(length(NonRedundantVals)/2)+NumElements];

    jack_f_180_0(1,iter) =  mean(NonRedundantVals(Back)) / mean(NonRedundantVals(For));
end

FreABfBo_std_f_180_0 = std(jack_f_180_0);

clear jack_f_180_0 For Back AnglesMinMinThres



%% Prepare outputStruct to save free segment and f_180_0 results
outputStruct(1).FreAfBo_mean_f_180_0 = FreAfBo_Single_f_180_0;
outputStruct(1).FreAfBo_std_f_180_0 = FreAfBo_std_f_180_0;
outputStruct(1).FreBfBo_mean_f_180_0 = FreBfBo_Single_f_180_0;
outputStruct(1).FreBfBo_std_f_180_0 = FreBfBo_std_f_180_0;
outputStruct(1).FreABfBo_mean_f_180_0 = FreABfBo_Single_f_180_0;
outputStruct(1).FreABfBo_std_f_180_0 = FreABfBo_std_f_180_0;
outputStruct(1).FreeSeg_AfBo = ffilt_FreeSeg_AfBo; % the index of Free segment after bound state segment, same as 'BoUnbo_CellTrackViterbiClass'
outputStruct(1).FreeSeg_BfBo = ffilt_FreeSeg_BfBo; % the index of Free segment before bound state segment
outputStruct(1).BoSegAfterFreeSeg_BfB = ffilt_BoSegAfterFreeSeg_BfB; % Consecutive bound segments after the free segments (those defined in FreeSeg_BfBo)
outputStruct(1).IdxLUT_FreeSeg_AfBo = ffilt_IdxLUT_FreeSeg_AfBo; % the same length as 'ffilt_FreeSeg_AfBo', but tell you the index of track index in 'BoUnbo_CellTrackViterbiClass'
outputStruct(1).IdxLUT_FreeSeg_BfBo = ffilt_IdxLUT_FreeSeg_BfBo; % the same length as 'ffilt_FreeSeg_BfBo', but tell you the index of track index in 'BoUnbo_CellTrackViterbiClass'
outputStruct(1).angleMatrix_FreAfBo = table(angleMatrix_FreAfBo(:,1),angleMatrix_FreAfBo(:,2),angleMatrix_FreAfBo(:,3),angleMatrix_FreAfBo(:,4),angleMatrix_FreAfBo(:,5),...
    'VariableNames',{'BoUnboTrackIndex','FirstJumpLength','SecJumpLength','Angle','CompleAngle'}); % table contains free segments after bound: the index of track index in 'BoUnbo_CellTrackViterbiClass', first segment length, second segment length, angles in rad, complementary angles in rad
outputStruct(1).angleMatrix_FreBfBo = table(angleMatrix_FreBfBo(:,1),angleMatrix_FreBfBo(:,2),angleMatrix_FreBfBo(:,3),angleMatrix_FreBfBo(:,4),angleMatrix_FreBfBo(:,5),...
    'VariableNames',{'BoUnboTrackIndex','FirstJumpLength','SecJumpLength','Angle','CompleAngle'}); % table contains free segments before bound: the index of track index in 'BoUnbo_CellTrackViterbiClass', first segment length, second segment length, angles in rad, complementary angles in rad

%% AUXILIARY FUNCTION %%

function ifBoUnboMix = CheckUnique(x)
    State = unique(x);
     
    if sum(State) == 3 % Bound/unbound segment mixed 
        ifBoUnboMix = true;
    else
        ifBoUnboMix = false;
    end
end

end

