function MSD = MSD_HMM_analyzer(CellTracks, CellTrackViterbiClass, MSD_time, FreeState)
%MSD_HMM_ANALYZER 
%   written by Anders Sejr Hansen (AndersSejrHansen@post.harvard.edu;
%   @Anders_S_Hansen; https://anderssejrhansen.wordpress.com)
%   License: GNU GPL v3
%   This function computes the time-averaged and ensemble-averaged MSD on
%   trajectories that have been HMM-classified, so it will only compute the
%   MSD on the free part. 

% store all displacements in a cell array; 
displacements = cell(1,length(MSD_time));

% the MSD should be a row vector
MSD = zeros(1,length(MSD_time));

% now go through each trajectory and consider only the free state:
for iter = 1:length(CellTracks)
           
    % the trajectory needs to have at least 1 displacement and at least 1
    % displacement between free states:
    if length(find(CellTrackViterbiClass{iter} == FreeState)) >= 1
        %%%% FIND THE RANGES OF INDICES IN WHICH THE PARTICLE IS FREE %%%%
        
        % find the jumps in the FreeState
        FreeIdx = find(CellTrackViterbiClass{iter} == FreeState);
        ranges = cell(0); %make a cell array to store the vectors
        
        %Now comes the complicated part: there can be multiple switches
        %between the free and bound states. So you need to figure out the
        %indices of these switches to ensure that you are only counting
        %adjacent Free parts.
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
                ranges{1} = FreeIdx(1:StateSwitches); % changed from FreeIdx(1):FreeIdx(StateSwitches);
                ranges{2} = FreeIdx(StateSwitches+1:end); % changed from FreeIdx(StateSwitches+1):FreeIdx(end)
            elseif length(StateSwitches) > 1
                %There are now at least three segments of free particles
                ranges{1} = FreeIdx(1:StateSwitches(1));
                for k=1:length(StateSwitches)-1
                    ranges{k+1} = FreeIdx(StateSwitches(k)+1:StateSwitches(k+1)); % changed from FreeIdx(StateSwitches(k)+1):FreeIdx(StateSwitches(k+1));
                end
                ranges{length(ranges)+1} = FreeIdx(StateSwitches(end)+1:end); % changed from FreeIdx(StateSwitches(end)+1):FreeIdx(end);
            end
        end
        
        %%%% COMPUTE THE MSD FROM DISPLACEMENTS IN THE FREE STATE %%%%
        for RangeIter = 1:length(ranges)
            % loop throigh all of the ranges of free displacements
            % CurrFreeIdx contain indices of the free displacements. So it
            % has length n, but need length n+1 of xy matrix:
            CurrFreeIdx = [ranges{RangeIter}; max(ranges{RangeIter})+1];
            curr_xy = CellTracks{iter}(CurrFreeIdx,:); % xy coordinates in microns      

            % need to do the calculations for all Delta t:
            AllDistancesMatrix = squareform(pdist(curr_xy, 'euclidean'));
            % This gives a matrix of distances:
            % [ 1,1     1,2     1,3     1,4     1,5 ]
            % [ 2,1     2,2     2,3     2,4     2,5 ]
            % [ 3,1     3,2     3,3     3,4     3,5 ]
            % [ 4,1     4,2     4,3     4,4     4,5 ]
            % [ 5,1     5,2     5,3     5,4     5,5 ]            
            % The diagonal of the matrix is useless
            % Need to get below the diagonal: use function diag:
            % e.g. diag(matrix, 1) gives: [2,1  3,2  4,3  5,4] 
            % e.g. diag(matrix, 4) gives: [5,1] 
            
            % find the min of trajectory length and timepoints:
            for DistIter = 1:min([length(MSD_time),  length(CurrFreeIdx)-1])
                % save all the displacements
                if isempty(displacements{1,DistIter})
                    displacements{1,DistIter} = diag(AllDistancesMatrix, DistIter);
                else
                    displacements{1,DistIter} = vertcat(displacements{1,DistIter}, diag(AllDistancesMatrix, DistIter));
                end
            end
        end
    end    
end


% You have now calculated all of the displacements between the free states.
% Now compute the MSD:
for iter = 1:length(MSD_time)
    curr_displacements = displacements{1,iter};
    MSD(1,iter) = mean(curr_displacements.^2);
end


end

