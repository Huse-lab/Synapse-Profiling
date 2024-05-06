function curvatureStats = calculateSimulatedCurvatureStats(z_data)

[mc, surfaceArea] = calculate_grid(z_data);

mc_v = mc(:);
mc_v = mc_v(mc_v~=0);

curvatureStats.mc_all = mc;
curvatureStats.SurfaceArea = surfaceArea;
curvatureStats.mean = mean(mc_v);
curvatureStats.range = range(mc_v);
curvatureStats.variance = var(mc_v);
curvatureStats.max_convexity = max(mc_v);
curvatureStats.mean_convexity = mean(mc_v(mc_v>0));
curvatureStats.total_convex_area = sum(mc_v>0)*(0.2^2); % 0.2 um/pixel or 0.04µm^2/pixel
curvatureStats.convex_area_map = mc > 0;


curvatureStats.max_concavity = min(mc_v);
curvatureStats.mean_concavity = mean(mc_v(mc_v<0));
curvatureStats.total_concave_area = sum(mc_v<0)*(0.2^2); % 0.2 um/pixel
curvatureStats.concave_area_map = mc < 0;

curvatureStats.PatternComplexity = calculateZernikeComplexity(mc,0.8);
curvatureStats.ZCoefficients = calculateZernikeCoefficients(mc,15);
curvatureStats.AdjustedZCoefficients = deRotateZernikeCoefficients(curvatureStats.ZCoefficients);

actin_map = curvatureStats.concave_area_map | z_data == 0;
curvatureStats.DegranMap = findEmptySlot(actin_map,5);
curvatureStats.DegranArea = sum(curvatureStats.DegranMap(:))*(0.2^2);
curvatureStats.TotalFreeArea = sum(~actin_map(:))*(0.2^2);

end


function [output_grid,surface_area] = calculate_grid(z_data)

    y_dim = size(z_data,1);
    x_dim = size(z_data,2);
    x1 = (1:x_dim)*(10000/x_dim);
    y1 = (1:y_dim)*(10000/y_dim);

    x1 = x1/1000;
    y1 = y1/1000;
    z = z_data/1000;
    [X,Y] = meshgrid(x1,y1);

    C = [];
    for j = 1:x_dim
        for i = 1:y_dim
            C = [C; X(i,j),Y(i,j),z(i,j)];
        end
    end

%     Z = griddata(x1,y1,z,X,Y);
    tri = delaunay(C(:,1),C(:,2));
    F.faces = tri;
    F.vertices = [C(:,1),C(:,2),C(:,3)];
    pc = GetCurvatures(F,0);
%     gc = pc(1,:).*pc(2,:);
    mc = -mean(pc);
    v1 = C(tri(:,2), :) - C(tri(:,1), :);
    v2 = C(tri(:,3), :) - C(tri(:,2), :);
    cp = 0.5*cross(v1,v2);
    surface_area = sum(sqrt(dot(cp, cp, 2)));



    output_grid = griddata(C(:,1),C(:,2),mc,X,Y);

end