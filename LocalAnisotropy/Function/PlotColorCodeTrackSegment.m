function PlotColorCodeTrackSegment(BoUnbo_CellTracks_idx,BoUnbo_CellTracks,BoUnbo_CellTrackViterbiClass,PlotWidth)
%PlotColorCodeTrackSegment Plot the turbo color coded segments of
%trajectory, where the start to end is shown as blue to red, and black
%segment represent the free segment, and colorful segment represent the
%bound segments.
%   Detailed explanation goes here
%   Zuhui Wang 2021/09/16


%% Plot segment colorcoded trajectory

TrackSegPos = BoUnbo_CellTracks{BoUnbo_CellTracks_idx};

free_segment = BoUnbo_CellTrackViterbiClass{BoUnbo_CellTracks_idx} == 2;

hold on
c = linspace(50,200,length(TrackSegPos(:,1)));
t = color_line(TrackSegPos(:,1),TrackSegPos(:,2),c,'LineWidth',1); %2
scatter(TrackSegPos(:,1),TrackSegPos(:,2),5,c,'filled'); %15
colormap(gca,'turbo')

title(['Track ID:' num2str(BoUnbo_CellTracks_idx)]);
% c = colorbar;
% c.Ticks = [1 50 100 150 200 250];
% c.TickLabels = {'start','','','','','end'};
line_color = t.CData; % 1 to end-1 row are effective color
line_color(free_segment,:) = repmat([1 1],sum(free_segment),1);
t.CData = line_color;
box on
ax = gca;
XLim = ax.XLim;
axis equal
ax.XLim = [mean(XLim)-PlotWidth/2, mean(XLim)+PlotWidth/2];
YLim = ax.YLim;
ax.YLim = [mean(YLim)-PlotWidth/2, mean(YLim)+PlotWidth/2];
hold off
end

