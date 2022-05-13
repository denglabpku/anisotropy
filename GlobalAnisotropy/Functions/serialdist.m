function [Distance] = serialdist(pair_position)
%serialdist calculate serial distance of paired position. 
%   Return distance between (x1,y1) & (x2,y2), (x2,y2) & (x3,y3), etc.
% pair_position:
%     [x_pos, y_pos]

D = pdist(pair_position);
Z = squareform(D);
Distance = diag(Z,1);
end

