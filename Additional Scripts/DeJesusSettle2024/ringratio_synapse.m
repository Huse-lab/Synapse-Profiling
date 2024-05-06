%% Actin radial profile from compiled_data
% 2024-01-03
% Analogue of "ring ratio" from discussion with AHS, stole function from
% analyze_curvaturedistancedata_v4.m

% NOTES FROM DISCUSSION WITH AHS:
%   1. Cut off the last thin bin, because when I draw ROIs based on actin,
%   I draw on the black space just outside the green signal, NOT on the
%   green signal itself -> explains the sudden drop off at the edges. After
%   checking different binnings of the distance-data structure, the outer
%   10 % of data is where the sudden dropoff happens. 2. Because of the
%   behavior of ratios, a "perfect ring ratio" is infinity; the scaling
%   behavior is exquisitely sensitive to small signals, which is
%   undesirable. Therefore, INVERSE RING RATIO (inside MFI / outside MFI)
%   is better and more stable. 3. Choosing the threshold of 0.7071 splits
%   the ratio into "inside" and "outside" based on keeping inside/outside
%   areas equal. Scale it according to the fraction of radial distances
%   taken. (ex. If inner 90% of radii are included, scale as 0.7071 x 0.9 =
%   0.6364.

%   Choose what fraction of the radii to include.
r_max = 0.9;
%   What should the ratio of the inner circle area over the outer annulus's be?
area_ratio = 2;

in_out_threshold = sqrt((area_ratio*(r_max)^2)/(1+area_ratio));
for idx=1:length(compiled_data)
    radial_distances = compiled_data(idx).Distance_NormalizedtoEdge;
    actin_signal = compiled_data(idx).CellStain;
    % Remove edge data
    within_edge = radial_distances < r_max;
    radial_distances = radial_distances(within_edge);
    actin_signal = actin_signal(within_edge);
    % IN / OUT
    ring_ratio = mean(actin_signal(...
        radial_distances <= in_out_threshold)) / ...
        mean(actin_signal(...
        radial_distances > in_out_threshold));
    compiled_data(idx).Ring_ratio_recip = ring_ratio;
end
