function curvatureData  = calculateCurvatureDistribution(MPStats,ROIdata)
%
%
%
% Alex Settle & Miguel de Jesus
% Memorial Sloan Kettering Cancer Center
% Morgan Huse Laboratory, Department of Immunology
% 2023


iMP=1;
[edgecoor,edgecoor_sph,mc,gc] = calculateCurvatures(MPStats(iMP));

theta = edgecoor_sph(:,1); phi = edgecoor_sph(:,2);
x = ROIdata(iMP).position(:,1); y = ROIdata(iMP).position(:,2); 
contact_map = inpolygon(theta,phi,x,y); 

edgecoor_sph_contact = edgecoor_sph(contact_map,:);
edgecoor_contact = edgecoor(contact_map,:);



centroid =[mean(edgecoor_sph_contact(:,1)),mean(edgecoor_sph_contact(:,2))];
distances = sqrt((edgecoor_sph_contact(:,1)'-centroid(1)).^2 + (edgecoor_sph_contact(:,2)'-centroid(2)).^2);

% 2. Execute computations per curvature value
roi_coords = ROIdata(iMP).position;
roi_coords(end+1,:) = roi_coords(1,:);
roi_shape = polyshape(roi_coords);
intersections = [];
for idx=1:length(edgecoor_sph_contact)
    endoflinex=centroid(1); endofliney=centroid(2); p=1;
    while isinterior(roi_shape,endoflinex,endofliney)
        endoflinex = centroid(1) + (edgecoor_sph_contact(idx,1)-centroid(1))*10^(p);
        endofliney = centroid(2) + (edgecoor_sph_contact(idx,2)-centroid(2))*10^(p);
        p=p+1;
    end
    % reshape roi_position data
    intersection_atcelledge = InterX([centroid(1) endoflinex;centroid(2) endofliney],[roi_coords(:,1)';roi_coords(:,2)']);
    if isempty(intersection_atcelledge)
        intersections(idx,1) = nan;
        intersections(idx,2) = nan;
    else
        intersections(idx,1) = intersection_atcelledge(1);
        intersections(idx,2) = intersection_atcelledge(2);
    end
    
end

% 3. Compute distances to cell edge
distances_tocelledge = sqrt((intersections(:,1)-centroid(1)).^2 + (intersections(:,2)-centroid(2)).^2);  
d_norm = distances'./distances_tocelledge;

curvatureData = struct();
curvatureData.edgecoor_contact = edgecoor_contact;
curvatureData.edgecoor_sph_contact = edgecoor_sph_contact;
curvatureData.MeanCurvature = mc(contact_map)';
curvatureData.GaussianCurvature = gc(contact_map)';
curvatureData.Distance_NormalizedtoEdge = d_norm;




end
