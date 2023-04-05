function newStats = scaleFeatures(featureStats,mc_i,scale)
%
% Function to apply a scaling factor to the protrusion or peak stats that
% are output by featurePicker.m. These will be in pixel values and this
% function will convert them to µm values in real space
%
% Alex Settle 
% Memorial Sloan Kettering Cancer Center
% Morgan Huse Laboratory, Department of Immunology
% 2023

newStats = featureStats;
N = size(mc_i,1);
scale_factor = scale/N; % µm/pixel
for k = 1:length(newStats)
    newStats(k).Area = newStats(k).Area*(scale_factor^2);
    newStats(k).Centroid = newStats(k).Centroid*scale_factor;
    newStats(k).BoundingBox = newStats(k).BoundingBox*scale_factor;
    newStats(k).MajorAxisLength = newStats(k).MajorAxisLength*scale_factor;
    newStats(k).MinorAxisLength = newStats(k).MinorAxisLength*scale_factor;
    newStats(k).ConvexHull = newStats(k).ConvexHull*scale_factor;
    newStats(k).ConvexArea = newStats(k).ConvexArea*(scale_factor^2);
    newStats(k).FilledArea = newStats(k).FilledArea*(scale_factor^2);
    newStats(k).Extrema = newStats(k).Extrema*scale_factor;
    newStats(k).EquivDiameter = newStats(k).EquivDiameter*scale_factor;
    newStats(k).Perimeter = newStats(k).Perimeter*scale_factor;
    newStats(k).PerimeterOld = newStats(k).PerimeterOld*scale_factor;
    newStats(k).WeightedCentroid = newStats(k).WeightedCentroid*scale_factor;
    newStats(k).MaxFeretDiameter = newStats(k).MaxFeretDiameter*scale_factor;
    newStats(k).MaxFeretCoordinates = newStats(k).MaxFeretCoordinates*scale_factor;
    newStats(k).MinFeretDiameter = newStats(k).MinFeretDiameter*scale_factor;
    newStats(k).MinFeretCoordinates = newStats(k).MinFeretCoordinates*scale_factor;

end




end