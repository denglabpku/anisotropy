function [filt_FreeSeg_AfBo,filt_FreeSeg_BfBo,filt_IdxLUT_FreeSeg_AfBo,filt_IdxLUT_FreeSeg_BfBo] = ...
    FindIdxofFreeSegmentInMix(BoUnbo_CellTracks,BoUnbo_CellTrackViterbiClass,FreSeg_size)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


% Find the free segment after or before a bound state segment
BoundState = 1;
FreeState = 2;
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
end

