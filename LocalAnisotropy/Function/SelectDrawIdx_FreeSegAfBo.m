function [FreAfBo_Single_f_180_0, Plot_CellTracks_idx,Plot_segment_idx] = SelectDrawIdx_FreeSegAfBo(sampleIter,FinalResults,DrawParams,BoUnbo_CellTracks)
%SelectDrawIdx_FreeSegAfBo Select the track index and corresponding
%segments indexes of one bound followed by free segments (bound -> free)
%   The function will return: FreAfBo_Single_f_180_0, f_180_0 of the
%   randomly selected segments; Plot_CellTracks_idx, index of
%   BoUnbo_CellTracks that containing bound->free segments;
%   Plot_segment_idx, the free segment index in each Plot_CellTracks_idx
%   selected BoUnbo_CellTracks.

    ContinueFreeSegments = DrawParams.ContinueFreeSegments ;
    PlotNumber = DrawParams.PlotNumber; % tptal number of visualized randomly selected trajectories
    MinMinJumpThres = DrawParams.MinMinJumpThres;
    MaxJump = DrawParams.MaxJump;   
%     FreeSeg_AfBo = FinalResults(4,sampleIter).FreeSeg_AfBo;
%     IdxLUT_FreeSeg_AfBo = FinalResults(4,sampleIter).IdxLUT_FreeSeg_AfBo;
    FreeSeg_AfBo = FinalResults(ContinueFreeSegments-1,sampleIter).FreeSeg_AfBo;
    IdxLUT_FreeSeg_AfBo = FinalResults(ContinueFreeSegments-1,sampleIter).IdxLUT_FreeSeg_AfBo;

    
    % Index of such FreeSeg_AfBo that free segment length are valid
    LengVal_FreeSeg_AfBo = [];
    
    % Calculate the jump length and discard tracks that not have valid
    % track length as in angle calcualtion step
    k = 1;
    kk = 1;
    for i = 1:length(IdxLUT_FreeSeg_AfBo)
        if length(FreeSeg_AfBo{i}) == ContinueFreeSegments
            trackStartIdx = FreeSeg_AfBo{i}(1);
            % trackEndIdx = FreeSeg_AfBo{i}(end);
            
            val_angle = 0;
            for j = 1:length(FreeSeg_AfBo{i})-1

                p1 = [BoUnbo_CellTracks{IdxLUT_FreeSeg_AfBo(i)}(trackStartIdx+j-1,1),BoUnbo_CellTracks{IdxLUT_FreeSeg_AfBo(i)}(trackStartIdx+j-1,2)];
                p2 = [BoUnbo_CellTracks{IdxLUT_FreeSeg_AfBo(i)}(trackStartIdx+j-1+1,1),BoUnbo_CellTracks{IdxLUT_FreeSeg_AfBo(i)}(trackStartIdx+j-1+1,2)];
                p3 = [BoUnbo_CellTracks{IdxLUT_FreeSeg_AfBo(i)}(trackStartIdx+j-1+2,1),BoUnbo_CellTracks{IdxLUT_FreeSeg_AfBo(i)}(trackStartIdx+j-1+2,2)];

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
                       angleMatrix_FreAfBo(k,:) = [IdxLUT_FreeSeg_AfBo(i), Distances, angle];
                       val_angle = val_angle+1;
                       k = k + 1;
                    end
                end
            end
            if val_angle == length(FreeSeg_AfBo{i})-1
                LengVal_FreeSeg_AfBo(kk) = i;
                kk = kk + 1;
            end
        end
    end
   
    Val_FreeSeg_AfBo = LengVal_FreeSeg_AfBo;
    
    % Generate the 16 random index of ValFreeSeg_AfBo
    RandVal_FreeSeg_AfBo = Val_FreeSeg_AfBo(randperm(length(Val_FreeSeg_AfBo),PlotNumber)); 

    Plot_CellTracks_idx = IdxLUT_FreeSeg_AfBo(RandVal_FreeSeg_AfBo); % Idx of to-be-plot track as BoUnbo_CellTrackViterbiClass

    Plot_segment_idx = FreeSeg_AfBo(RandVal_FreeSeg_AfBo); % Idx of free-segments-to-be-plot with track idx as Plot_CellTracks_idx
    
    %% Calculate f_180_0 of free segments after bound segment

    nBins = 1*36; %so 10 degrees per bin
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

end

