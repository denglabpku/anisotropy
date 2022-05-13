function [FreBfBo_Single_f_180_0,Plot_CellTracks_idx,Plot_segment_idx] = SelectDrawIdx_FreeSegBfBo(sampleIter,FinalResults,DrawParams,BoUnbo_CellTracks)
%SelectDrawIdx_FreeSegBfBo Similar to SelectDrawIdx_FreeSegAfBo, but
%randomly selecting free -> bound segments


    ContinueFreeSegments = DrawParams.ContinueFreeSegments ;
    PlotNumber = DrawParams.PlotNumber; % tptal number of visualized randomly selected trajectories
    MinMinJumpThres = DrawParams.MinMinJumpThres;
    MaxJump = DrawParams.MaxJump;   
%     FreeSeg_BfBo = FinalResults(4,sampleIter).FreeSeg_BfBo;
%     IdxLUT_FreeSeg_BfBo = FinalResults(4,sampleIter).IdxLUT_FreeSeg_BfBo;
    FreeSeg_BfBo = FinalResults(ContinueFreeSegments-1,sampleIter).FreeSeg_BfBo;
    IdxLUT_FreeSeg_BfBo = FinalResults(ContinueFreeSegments-1,sampleIter).IdxLUT_FreeSeg_BfBo;
    
    % Index of such FreeSeg_AfBo that free segment length are valid
    LengVal_FreeSeg_BfBo = [];
    
    % Calculate the jump length and discard tracks that not have valid
    % track length as in angle calcualtion step
    k = 1;
    kk = 1;
    for i = 1:length(IdxLUT_FreeSeg_BfBo)
        if length(FreeSeg_BfBo{i}) == ContinueFreeSegments
            trackStartIdx = FreeSeg_BfBo{i}(1);
            % trackEndIdx = FreeSeg_AfBo{i}(end);
            
            val_angle = 0;
            for j = 1:length(FreeSeg_BfBo{i})-1

                p1 = [BoUnbo_CellTracks{IdxLUT_FreeSeg_BfBo(i)}(trackStartIdx+j-1,1),BoUnbo_CellTracks{IdxLUT_FreeSeg_BfBo(i)}(trackStartIdx+j-1,2)];
                p2 = [BoUnbo_CellTracks{IdxLUT_FreeSeg_BfBo(i)}(trackStartIdx+j-1+1,1),BoUnbo_CellTracks{IdxLUT_FreeSeg_BfBo(i)}(trackStartIdx+j-1+1,2)];
                p3 = [BoUnbo_CellTracks{IdxLUT_FreeSeg_BfBo(i)}(trackStartIdx+j-1+2,1),BoUnbo_CellTracks{IdxLUT_FreeSeg_BfBo(i)}(trackStartIdx+j-1+2,2)];

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
                       angleMatrix_FreBfBo(k,:) = [IdxLUT_FreeSeg_BfBo(i), Distances, angle];
                       val_angle = val_angle+1;
                       k = k+1;
                    end
                end
            end
            if val_angle == length(FreeSeg_BfBo{i})-1
                LengVal_FreeSeg_BfBo(kk) = i;
                kk = kk + 1;
            end
        end
    end
   
    Val_FreeSeg_BfBo = LengVal_FreeSeg_BfBo;
    
    % Generate the 16 random index of ValFreeSeg_AfBo
    RandVal_FreeSeg_BfBo = Val_FreeSeg_BfBo(randperm(length(Val_FreeSeg_BfBo),PlotNumber)); 

    Plot_CellTracks_idx = IdxLUT_FreeSeg_BfBo(RandVal_FreeSeg_BfBo); % Idx of to-be-plot track as BoUnbo_CellTrackViterbiClass

    Plot_segment_idx = FreeSeg_BfBo(RandVal_FreeSeg_BfBo); % Idx of free-segments-to-be-plot with track idx as Plot_CellTracks_idx
    
    %% Calculate f_180_0 of free segments before bound segment
    nBins = 1*36; %so 10 degrees per bin
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

end

