function [ trackedPar_merged_QC ] = RemoveAmbigiousTracks( trackedPar_merged, ClosestDist )
%REMOVEAMBIGIOUSCONNECTIONS 
%   written by Anders Sejr Hansen (AndersSejrHansen@post.harvard.edu;
%   @Anders_S_Hansen; https://anderssejrhansen.wordpress.com/)
%   License: GNU GPL v3

%   Go through all trajectories and remove ambigious tracks

%   STRATEGY
%   make a list of all locs in each frame
%   go through each trajectory
%   calculate distance to each other loc
%   If two locs closer than ClosestDist, then stop both trajectories there
%   and discard the rest
%   Repeat for all trajectories

%%%%%%%%%%%% Make a list for all localizations: %%%%%%%%%%%%
final_frame = max(trackedPar_merged(1,end).Frame);
AllLocs = cell(1,final_frame+100);
%   All Locs has 3 columns:
%   Column 1: x
%   Column 2: y
%   Column 3: trajID (element in trackedPar)

for iter=1:length(trackedPar_merged)
    curr_xy = trackedPar_merged(1,iter).xy;
    curr_frames = trackedPar_merged(1,iter).Frame;
    for sub_iter = 1:length(curr_frames)
        AllLocs{1,curr_frames(sub_iter)} = vertcat(AllLocs{1,curr_frames(sub_iter)}, [curr_xy(sub_iter,:), iter]);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%% Loop through and remove ambigious tracks %%%%%%%%%%%%
trackedPar_merged_QC = struct();
QC_counter = 1;
% loop through trackedPar:
for iter=1:length(trackedPar_merged) %iter through every tracks
    curr_xy = trackedPar_merged(1,iter).xy;
    curr_frames = trackedPar_merged(1,iter).Frame;
    % if there is only 1 localization, you can just save it:
    if length(curr_frames) == 1
        trackedPar_merged_QC(1,QC_counter).xy = trackedPar_merged(1,iter).xy;
        trackedPar_merged_QC(1,QC_counter).Frame = trackedPar_merged(1,iter).Frame;
        trackedPar_merged_QC(1,QC_counter).TimeStamp = trackedPar_merged(1,iter).TimeStamp;
        QC_counter = QC_counter + 1; % increment the counter;
    else
        % OK there was at least one displacement / tracking connection. So
        % check to see if it was OK or not:
        MaxFrame = length(curr_frames); % assume you make it to the end
        for sub_iter = 1:length(curr_frames) %sub_iter thru each frames in the iter-th track
            ThisFrame = curr_frames(sub_iter);
            ThisTrajID = iter;
            ThisXY = curr_xy(sub_iter,:);
            LocsThisFrame = AllLocs{1,ThisFrame};
            % find all locs other than the one from this trajectory:
            LocsToConsider = find(LocsThisFrame(:,3) ~= ThisTrajID);
            if ~isempty(LocsToConsider)
                % calculate distance to other locs:
                curr_dists = zeros(1,length(LocsToConsider));
                for subsub_iter = 1:length(LocsToConsider) %subsub_iter thru each other locs in the sub_iter-th frame
                    curr_dists(1,subsub_iter) = sqrt( (ThisXY(1,1)-LocsThisFrame(LocsToConsider(subsub_iter),1))^2  +   (ThisXY(1,2)-LocsThisFrame(LocsToConsider(subsub_iter),2))^2   );
                end
                
                %Check to see if you got too close to another particle
                if min(curr_dists) < ClosestDist
                    % Oops: you got too close to another particle; cannot
                    % reliably do the tracking. Need to stop tracking and
                    % throw away the rest of the trajectory
                    MaxFrame = max([1, sub_iter-1]);
                    break; 
                end
            end
        end
        % OK, now you have determined which was the bigest acceptable frame
        % to continue the trajectory, MaxFrame. So save the QC'ed
        % trajectory:
        trackedPar_merged_QC(1,QC_counter).xy = trackedPar_merged(1,iter).xy(1:MaxFrame,:);
        trackedPar_merged_QC(1,QC_counter).Frame = trackedPar_merged(1,iter).Frame(1:MaxFrame);
        trackedPar_merged_QC(1,QC_counter).TimeStamp = trackedPar_merged(1,iter).TimeStamp(1:MaxFrame);
        QC_counter = QC_counter + 1; % increment the counter;
        
    end    
end