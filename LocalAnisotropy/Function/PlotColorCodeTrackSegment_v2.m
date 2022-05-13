function PlotColorCodeTrackSegment_v2(Plot_CellTracks_idx,Plot_segment_idx,BoUnbo_CellTracks,DrawParams)
%PlotColorCodeTrackSegment_v2 
%%  Description
% Plot the turbo color coded segments of
% trajectory, where the start to end is shown as blue to red, and black
% segment represent the free segment, and colorful segment represent the
% bound segments.The red dashline circle centered around the last locs of
% bound segments with the defined radius (units: um).
%   Dependency: Colored line or scatter plot (v1.0.0.0) from Matlab Add-Ons
    % Pekka Kumpulainen (2021). Colored line or scatter plot
    % (https://www.mathworks.com/matlabcentral/fileexchange/19476-colored-line-or-scatter-plot),
    % MATLAB Central File Exchange. Retrieved October 14, 2021.
%   Zuhui Wang 2021/09/16


%% Plot segment colorcoded trajectory
PlotWidth = DrawParams.PlotWidth; % unit: um
PlotMode = DrawParams.PlotMode ; % tptal number of visualized randomly selected trajectories
CirRadius = DrawParams.CircleRadius; % units: um

if PlotMode == 1 % if draw segment: bound -> free
    TrackSegPos = BoUnbo_CellTracks{Plot_CellTracks_idx};
    bound_segment = Plot_segment_idx(1)-1;
    free_segment = Plot_segment_idx;
    
    draw_segment = [bound_segment free_segment];
    draw_pos = [draw_segment draw_segment(end)+1];
    
    hold on
    c_freeline = ones(1,length(free_segment)+1);
    c_line = [50,c_freeline]; % start from blue bound to black free segments
    c_pos = linspace(50,175,length(draw_pos));
    color_line(TrackSegPos(draw_pos,1),TrackSegPos(draw_pos,2),c_line,'LineWidth',1); %2
    scatter(TrackSegPos(draw_pos,1),TrackSegPos(draw_pos,2),5,c_pos,'filled'); %15
    viscircles(TrackSegPos(draw_pos(2),:),CirRadius,'LineStyle','-','LineWidth',1);
    colormap(gca,'turbo')
    title(['Track ID:' num2str(Plot_CellTracks_idx)]);
    box on
    ax = gca;
    XLim = ax.XLim;
    axis equal
    ax.XLim = [mean(XLim)-PlotWidth/2, mean(XLim)+PlotWidth/2];
    YLim = ax.YLim;
    ax.YLim = [mean(YLim)-PlotWidth/2, mean(YLim)+PlotWidth/2];
    hold off
    
    
elseif PlotMode == 2 % if draw free segments before bound
    TrackSegPos = BoUnbo_CellTracks{Plot_CellTracks_idx};
    free_segment = Plot_segment_idx;
    bound_segment = Plot_segment_idx(end)+1;
    
    draw_segment = [free_segment bound_segment];
    draw_pos = [draw_segment draw_segment(end)+1];
    
    hold on
    c_freeline = ones(1,length(free_segment));
    c_line = [c_freeline,154,154]; % start from black free segments to red bound segments
    c_pos = linspace(50,175,length(draw_pos));
    color_line(TrackSegPos(draw_pos,1),TrackSegPos(draw_pos,2),c_line,'LineWidth',1); %2
    scatter(TrackSegPos(draw_pos,1),TrackSegPos(draw_pos,2),5,c_pos,'filled'); %15
    viscircles(TrackSegPos(draw_pos(end-1),:),CirRadius,'LineStyle','-','LineWidth',1);
    colormap(gca,'turbo')
    title(['Track ID:' num2str(Plot_CellTracks_idx)]);
    box on
    ax = gca;
    XLim = ax.XLim;
    axis equal
    ax.XLim = [mean(XLim)-PlotWidth/2, mean(XLim)+PlotWidth/2];
    YLim = ax.YLim;
    ax.YLim = [mean(YLim)-PlotWidth/2, mean(YLim)+PlotWidth/2];
    hold off
end


end

