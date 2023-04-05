function [protrusionStats,peakStats,labeledProtrusions,labeledPeaks] = featurePicker(mc_grid)
% featurePicker.m -function for identifying features within
% 50x50 mean-curvature synapse grids. 
% 1/5/2022 start
% Alex Settle + Miguel de Jesus
% Memorial Sloan Kettering Cancer Center
% Morgan Huse Laboratory, Department of Immunology
% 2023
%
% Input arguments
%   - mc_grid - NxN mean curvature matrix
%
% Outputs
%   - protrusionStats - struct containing list of protrustions
%          (negative curvature valleys) with relevant statistics for each
%   - peakStats - struct containing list of peaks
%          (positive curvature ) with relevant statistics for each
%
%   - labeledProtrusions - a 100x100 grid with color-coded protrusions
%   labeled 
%   - labeledPeaks - same for peaks
%


%interpolate grid for smoother running
N = size(mc_grid,1);
[X,Y] = meshgrid(1:N);
[Xq,Yq] = meshgrid(1:.5:N);
%gaussian filter
mcq = interp2(X,Y,mc_grid,Xq,Yq);
mcfilt = imgaussfilt(mcq,3);

water = watershed(mcfilt);
inv_water = watershed(-mcfilt);

%Get local maxima and minima

neg_threshold = 0;
neg_mask = mcfilt < neg_threshold;

pos_threshold = 0.001;
pos_mask = mcfilt > pos_threshold;

labeledProtrusions = water;
labeledProtrusions(neg_mask == 0) = 0;

labeledPeaks = inv_water;
labeledPeaks(pos_mask == 0) = 0;

protrusionStats = regionprops(labeledProtrusions,mcq,'all');
protrusionStats = protrusionStats([protrusionStats.Area] > 5);

peakStats = regionprops(labeledPeaks,mcq,'all');
peakStats = peakStats([peakStats.Area] > 5);

return

end

