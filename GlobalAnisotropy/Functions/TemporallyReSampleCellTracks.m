function [ new_CellTracks, new_CellTrackViterbiClass ] = TemporallyReSampleCellTracks(CellTracks, CellTrackViterbiClass, TimePointIndex)
%TEMPORALLYRESAMPLETRACKEDPAR get other delta-time points
%   written by Anders Sejr Hansen (AndersSejrHansen@post.harvard.edu;
%   @Anders_S_Hansen; https://anderssejrhansen.wordpress.com)
%   License: GNU GPL v3
%   Take as input HMM-classified CellTracks

% if TimePointsIndex = 2, every other timepoint is used: e.g. 1,3,5,7.
% If TimePointsIndex = 4, every fourth is used: e.g. 1, 5, 9, 13

%%%%%%%%%%% Loop through and remove ambigious tracks %%%%%%%%%%%%
TempStructure = struct();
NextTrajIter = 1;

for iter = 1:length(CellTracks)
    % check to see if trajectory is long enough:
    % This is a bit complicated, because the trajectories may contain gaps
    % so the number of localizations in of themselves do not tell the whole
    % story. 
    curr_Frames = 1:size(CellTracks{iter},1);
    if length(curr_Frames) > TimePointIndex
        % OK, there is at least one displacement to sample:
        % loop through as follows:
        %   start with the first frame n, then find n+TimePointIndex,
        %   n+2TimePointIndex, n+3TimePointIndex, etc.
        %   Then proceed with the 2nd frame n+1, then find
        %   n+1+TimePointIndex, n+1+2*TimePointIndex, etc.
        % 
        %   But be careful to avoid overcounting. E.g. if TimePointIndex=2,
        %   then starting at frame 1 and frame 3 both will have many of the
        %   same displacements, e.g. 3-5 will be contained in both. So set 
        %   MaxStartIndex to the first frame + TimePointIndex - 1;
        MaxStartIndex = min([max(curr_Frames) - TimePointIndex, min(curr_Frames) + TimePointIndex - 1]);
        
        % now loop through:
        for IndexIter = 1:MaxStartIndex
            % since there are no gaps, it's very simple:
            chosen_Frames = IndexIter:TimePointIndex:length(curr_Frames);
            
            % OK, you finished compiling a trajectory. So now save it:
            temp_CellTracks_matrix = CellTracks{iter};
            TempStructure(1,1).CellTracks{1,NextTrajIter} = temp_CellTracks_matrix(chosen_Frames,:);
            
            % Now the HMM classification is a bit complicated: if you have
            % n localizations, you have n-1 HMM-classifications of jumps.
            % But if you take a long jump, the classification of all the
            % intermediate jumps may be a combination of FREE and BOUND.
            % Here we want to be conservative: we really do not want bound
            % molecules contaminating the free fraction: so set a very stringest threshold:
            % uncless all molecules are free, set it to bound:
            curr_Viterbi = CellTrackViterbiClass{iter};
            new_Viterbi = ones(length(chosen_Frames)-1,1);
            for HMM_iter = 1: (length(chosen_Frames)-1)
                first_idx = chosen_Frames(HMM_iter);
                last_idx = chosen_Frames(HMM_iter+1)-1; % subtract 1 since Viterbi is one element shorter than CellTracks
                curr_class = curr_Viterbi(first_idx:last_idx);
                if min(curr_class) == 2
                    new_Viterbi(HMM_iter,1) = 2;
                else
                    new_Viterbi(HMM_iter,1) = 1;
                end
            end
            % also save the HMM classification: 
            TempStructure(1,1).CellTrackViterbiClass{1,NextTrajIter} = new_Viterbi;
            % update the counter:
            NextTrajIter = NextTrajIter + 1;
        end
    end
end    

% OK, you are done: save as cell arrays instead of structured arrays:
new_CellTracks = TempStructure(1,1).CellTracks;
new_CellTrackViterbiClass = TempStructure(1,1).CellTrackViterbiClass;
end

