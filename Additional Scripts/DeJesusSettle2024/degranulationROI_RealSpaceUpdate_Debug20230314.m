% ROI selection script for degranulation (2022-11-17)

% Pseudo-algorithm:
% * Based on selectROI_v5.m
% Start with MPRenders, ask user for which frames of interest to process
% Plot R & cell channel maps for selecting ROI
% Three rounds of selection: 
%     1. ROI for the degranulation spot (cell channel)
%     - Click a spot of interest
%     - Create viscircles of a certain radius
%     - User confirm
%     - Store information 
%     2. ROI for non-degranulation spots (cell channel)
%     - Repeat #1 process
%     - Store in the same spot as the degran spot array index=2:end
%     3. ROI for the entire synapse (cell OR R channel)
%     - Draw freehand
%     - User confirm
%     - Store information as its own separate field
% Store information into structure for analysis
% * Based on analyze_curvaturedistancedata_v4.m
% Collect coordinates, mean curvatures, & Gauss curvatures in both polygons

% 2022-11-22: CODE STILL FAILS FOR REPEATING CONTROL POINTS because of
% loop-break error

% 2022-11-22: MH notes that can use MC map to select the inner-concave area
% as the "whole-synapse" control. Something to re-code.

% 2022-01-10-12: Re-structured the entire code to return to mapping based
% on displacement/shape, picking dots, & streamlining curvature collection
% based on REAL SPACE coordinates. (See final function)

%% Initialization

path = uigetdir(pwd,'Select folder containing MPRenders.');
file_list = dir(fullfile(path, '*MPRender*'));

warning('off','all')
shellsize = linspace(.1,1,10);

degranulation_topography = struct('Label',{},'Timestamp',[],...
'Spots_Degran',[],'Synapse_Coords_Sph',[],...
'MC_Degran',{},'MC_Synapse',{},...
'GC_Degran',{},'GC_Synapse',{});

%% Loop through MPStats files

% Open MPRender_pool IMG147 and ROI_pool IMG147 with the aggregated time
% points already

% For timepoints & selection criteria, refer to slide 1 of IM147_DegranulationEvents_Compiled.pptx
 
for idx=1:length(MPStats)
    
    disp(strcat(MPStats(idx).FileName," timepoint ",num2str(MPStats(idx).Timestamp)))
    [edgecoor,centroidcoor,edgecoor_sph,TRI_nosph] = check_processing_status(MPStats(idx));
    
    % Calculate curvatures
        theta = edgecoor_sph(:,1); phi = edgecoor_sph(:,2);
        ThetaPhiR = edgecoor_sph;         
        [ThetaPhiRsmooth] = Smooth_spherical_mesh(edgecoor_sph,sqrt(1^2/pi),'arc',1,'mean');
        [X,Y,Z] = sph2cart(ThetaPhiRsmooth(:,1),ThetaPhiRsmooth(:,2),ThetaPhiRsmooth(:,3));    
        FV.faces    = MPStats(idx).TRI_Connectivity_SPHbound;
        FV.vertices = [X,Y,Z];
        pc = GetCurvatures(FV,0);
        gc = pc(1,:).*pc(2,:);
        mc = mean(pc);
        
    % Skip synapse selection, use my pre-selected ROIs from Alex's code.            
        roi_synapse = roi_position(idx).position;
        
    % Loop within a synapse selection until degranulation points are OK.
        done_degran = false;
        while done_degran == false
            
        % Plot the 2D surface
            fprintf('Click on the degranulation spot/s and then press ENTER.\n');
            figure('Name',MPStats(idx).FileName,'position',[0,300,2000,700])

        % Intensity-based coloring
            stain_I = MPStats(idx).stain_int;
            contact_map = inpolygon(theta,phi,roi_synapse(:,1),roi_synapse(:,2));
            stain_I_withinROI = stain_I.*contact_map;
            ph_synapse = patch('faces',TRI_nosph,'vertices',[ThetaPhiR(:,[1 2]) stain_I_withinROI],'facevertexcdata',stain_I_withinROI,'facecolor','flat','edgecolor','none');
            caxis([prctile(stain_I,1) prctile(stain_I,99)]) 
            Topograph_Aesthetics
            % Define a custom green colormap for the intensity           
                stainmap = ones(256,3)*.05;
                stainmap(:,2) = linspace(.005,.995,256);
                colormap(stainmap)                        

        %Prompt user to select ROI for degranulation spot
            hold on
            [roi_degran_x,roi_degran_y] = getpts;
            roi_degran = {};
            x_bound = xlim;
            radius = (x_bound(2)-x_bound(1))/60;
                for jdx=1:length(roi_degran_x)
                    drawcircle('Center',[roi_degran_x(jdx),roi_degran_y(jdx)],'Radius',radius);
                end
        % Make sure the selections are satisfactory
            accept = questdlg('Are you satisfied with these degranulation spots?',MPStats(idx).FileName,'Yes','No, Try Again','Yes');
            if strcmp(accept,'Yes')
                done_degran = true;
                close
            else
                close               
            end
        end
       
% Collect curvature data within desired shell-size parameters, cleaning up output. 

    % Copy data into output structure.
        degranulation_topography(idx).Label = MPStats(idx).FileName;
        degranulation_topography(idx).Timestamp = MPStats(idx).Timestamp;
        degranulation_topography(idx).Spots_Degran = [roi_degran_x roi_degran_y];
        degranulation_topography(idx).Synapse_Coords_Sph = roi_synapse;
    % Function to calculate the coordinates within real-distance shells based on angular space selections. 
        curvatures_withinShell = GetCurvatures_withinShells(MPStats(idx),degranulation_topography(idx),shellsize,mc,gc);
    % Final line to copy over data into struct roi_degran_position...
    % Degranulation spots
        degranulation_topography(idx).MC_Degran = curvatures_withinShell.MC_Degran;
        degranulation_topography(idx).GC_Degran = curvatures_withinShell.GC_Degran;
    % Entire synapse
        degranulation_topography(idx).MC_Synapse = curvatures_withinShell.MC_Synapse;
        degranulation_topography(idx).GC_Synapse = curvatures_withinShell.GC_Synapse;

end

%% Minimal checking based on mis-alignment in Daan's function: with AHS 20230316

% Pre-load MPRender_pool & Degranulation Topography (just for the
% pre-selected degranulation spots)

shellsize = linspace(.1,1,10);

for idx=1:15
    
    % Initialize Daan curvature calcs
        [edgecoor,centroidcoor,edgecoor_sph,TRI_nosph] = check_processing_status(MPStats(idx));
    % Calculate curvatures
        theta = edgecoor_sph(:,1); phi = edgecoor_sph(:,2);
        ThetaPhiR = edgecoor_sph;         
        [ThetaPhiRsmooth] = Smooth_spherical_mesh(edgecoor_sph,sqrt(1^2/pi),'arc',1,'mean');
        [X,Y,Z] = sph2cart(ThetaPhiRsmooth(:,1),ThetaPhiRsmooth(:,2),ThetaPhiRsmooth(:,3));    
        FV.faces    = MPStats(idx).TRI_Connectivity_SPHbound;
        FV.vertices = [X,Y,Z];
        pc = GetCurvatures(FV,0);
        gc = pc(1,:).*pc(2,:);
        mc = mean(pc);
        
    % Select ROI based on curvature map without rim. Comment out if we
    % don't want that.
        degranulation_topography(idx).Synapse_Coords_Sph = roi_position(idx).position;
        
%     % Loop within a time point until the synapse selection is OK.
%         done_synapse = false;
%         while done_synapse==false
% 
%         % Plot the 2D surface
%             figure('Name',MPStats(idx).FileName,'position',[0,300,2000,700])
%             Topograph_Aesthetics
%             % Switch between the three lines below depending on what to plot
%             ph_synapse = patch('faces',TRI_nosph,'vertices',[ThetaPhiR(:,[1 2]) mc'],'facevertexcdata',mc','facecolor','flat','edgecolor','none');
%             colormap(turquoisebrown)
%             caxis([prctile(mc,1) prctile(mc,99)]) 
% 
%         % Prompt user to select ROI
%             roi_synapse = drawfreehand('FaceAlpha',0.1,'Color',[1 1 0]);        
%             accept = questdlg('Are you satisfied with this ROI selection?', MPStats(idx).FileName,'Yes','No, Try Again','Yes');
%             if strcmp(accept,'Yes')
%                 done_synapse = true;
%                 degranulation_topography(idx).Synapse_Coords_Sph = roi_synapse.Position;
%             else
%                 close               
%             end           
%         end

    % Function to calculate the coordinates within real-distance shells based on angular space selections. 
        curvatures_withinShell = GetCurvatures_withinShells(MPStats(idx),degranulation_topography(idx),shellsize,mc,gc);
    % Final line to copy over data into struct roi_degran_position...
    % Degranulation spots
        degranulation_topography(idx).MC_Degran = curvatures_withinShell.MC_Degran;
        degranulation_topography(idx).GC_Degran = curvatures_withinShell.GC_Degran;
    % Entire synapse
        degranulation_topography(idx).MC_Synapse = curvatures_withinShell.MC_Synapse;
        degranulation_topography(idx).GC_Synapse = curvatures_withinShell.GC_Synapse;

end

%% Calculate radial distance of degranulation events

list = string(shellsize);
[shellsize_choice,~] = listdlg('Promptstring',{'Select a spot size in microns:'},'ListString',list,'SelectionMode','single');
warning('off','all');

degranulation_radialized = struct('coord_spot',[],'coord_synapse',[],'MC_degran',[],'GC_degran',[]);
rowdata = degranulation_radialized;
% Copy over data with one row for individual degranulation spots.
for idx=1:length(degranulation_topography)
    for jdx=1:size(degranulation_topography(idx).Spots_Degran,1)
        rowdata.coord_spot = degranulation_topography(idx).Spots_Degran(jdx,:);
        rowdata.coord_synapse = degranulation_topography(idx).Synapse_Coords_Sph;
        % Pick radius-um local spot MC to average for the y-axis.
        rowdata.MC_degran = degranulation_topography(idx).MC_Degran{jdx,shellsize_choice};
        rowdata.GC_degran = degranulation_topography(idx).GC_Degran{jdx,shellsize_choice};
        if idx==1
            degranulation_radialized = rowdata;
        elseif idx>1
            degranulation_radialized = [degranulation_radialized rowdata];
        end
    end
end

for idx=1:length(degranulation_radialized)
    
    centroid = [mean(degranulation_radialized(idx).coord_synapse(:,1)),mean(degranulation_radialized(idx).coord_synapse(:,2))];
    spot = degranulation_radialized(idx).coord_spot;
    d_spot2centroid = norm(centroid-spot);
    
    % Find distance to cell edge crossing the spot.
   
    roi_coords = degranulation_radialized(idx).coord_synapse;
    roi_coords(end+1,:) = roi_coords(1,:);
    roi_shape = polyshape(roi_coords);
    intersections = [];
        endoflinex=centroid(1); endofliney=centroid(2); p=1;
        while isinterior(roi_shape,endoflinex,endofliney)
            endoflinex = centroid(1) + (spot(1)-centroid(1))*10^(p);
            endofliney = centroid(2) + (spot(2)-centroid(2))*10^(p);
            p=p+1;
        end
        % reshape roi_position data
        intersection_atcelledge = InterX([centroid(1) endoflinex;centroid(2) endofliney],[roi_coords(:,1)';roi_coords(:,2)']);
        if isempty(intersection_atcelledge)
            intersections(jdx,1) = nan;
            intersections(jdx,2) = nan;
        else
            intersections(jdx,1) = intersection_atcelledge(1);
            intersections(jdx,2) = intersection_atcelledge(2);
        end
    distances_tocelledge = sqrt((intersections(:,1)-centroid(1)).^2 + (intersections(:,2)-centroid(2)).^2);  
    d_norm = d_spot2centroid./distances_tocelledge;
    degranulation_radialized(idx).Normalized_D = d_norm;

end

%% Plot: Collect MC & GC and overlay degranulations

% If compiled_data is open
x_d = [compiled_data.Distance_NormalizedtoEdge];
y_mc = [compiled_data.MeanCurvature];
y_gc = [compiled_data.GaussianCurvature];

% Sliding window to collect y_mc, y_gc data & SDs around +/-0.05 x_d
    x_transform = []; y_transform = []; y_SD=[];
    for idx=1:length(x_d)
        window = [x_d(idx)-.05 x_d(idx)+.05];
        x_within = x_d(x_d>=window(1) & x_d<=window(2));
        y_within = y_mc(x_d>=window(1) & x_d<=window(2));
        % Take the sliding window data averages
        x_transform(idx) = mean(x_within);
        y_transform(idx) = mean(y_within);
        y_SD(idx) = std(y_within);
    end
    % Re-sort everybody
        [x_transform,sortidx] = sort(x_transform,'ascend');
        y_transform = y_transform(sortidx);
        y_SD = y_SD(sortidx);
            
    % Plot MC
%     subplot(2,1,1),hold on
    figure, hold on
    errorbar(x_transform,y_transform,y_SD,'Color',[.75 .75 .75],'CapSize',0);
    scatter(x_transform,y_transform,10,[0 0 0],'.');
    % Overlay degranulation spots
    spot_x = [degranulation_radialized.Normalized_D];
    for idx=1:length(degranulation_radialized)
        spot_y(idx) = mean(degranulation_radialized(idx).MC_degran);
    end
    scatter(spot_x,spot_y,100,[1 0 0],'x','LineWidth',2)
    ylim([-1 1]),daspect([1 4 1])
    xlabel('Radial distance'),ylabel('Mean curvature')
    
%% Plotting empirical CDF overlays

% Initialize y_mc from previous section
y_mc = [compiled_data.MeanCurvature];
y_spot = [degranulation_radialized.MC_degran];

figure, hold on
cdf_synapse = cdfplot(y_mc);
x_all = cdf_synapse.XData;
y_all = cdf_synapse.YData;
cdf_spot = cdfplot(y_spot);
x_spot = cdf_spot.XData;
y_spot = cdf_spot.YData;
xlim([-1 1]), pbaspect([2 1 1])
figure, hold on
plot(x_all,y_all,'LineWidth',5,'Color',[.4 .4 .4])
plot(x_spot,y_spot,'LineWidth',5,'Color',[1 0 0])
xlim([-1 1]), pbaspect([2 1 1]), xticks([-1:.5:1])

[h,p,D] = kstest2(y_mc,y_spot)
    
%% Functions

% Select particle for plot
function iMP = Select_particle_for_plot(MPStats,iMP)
    NMPs = length(MPStats);
    if NMPs > 1
        iMP = str2double(inputdlg('Which timepoint do you want to visualize?',...
            'Select particle',[1 50],{num2str(iMP)}));
        if iMP > NMPs
            error(['Particle ' num2str(iMP) ' doesn''t exist']);
        end
    else
        iMP = 1;
    end
end

% Set 2D plot options
function[] = Topograph_Aesthetics()

    axis off, axis equal, view(0,90);
    xlabel('Theta')
    ylabel('Phi')
    grid off
    set(gca,'Clipping','off')
    add_colorscalebar()
    colormap(turquoisebrown)

end

% Add colorbar to 2D plot, while keeping the same axis position
function[] = add_colorscalebar()

    ch = colorbar;
    chpos = get(ch,'Position');
    original_axis_pos  = get(gca,'Position');
    original_axis_pos(2) = (1-original_axis_pos(4))/2;
    set(ch,'Position',[chpos(1) 0.5-original_axis_pos(4)/4 chpos(3) original_axis_pos(4)/2]);
    set(gca,'Position',original_axis_pos);

end

function [xo yo R] = circle_fit(x,y)
% A function to find the best circle fit (radius and center location) to
% given x,y pairs
% 
% Val Schmidt
% Center for Coastal and Ocean Mapping
% University of New Hampshire
% 2012
% 
% Arguments:
% x:         x coordinates
% y:         y coordinates
%
% Output:
% xo:        circle x coordinate center
% yo:        circle y coordinate center
% R:         circle radius

x = x(:);
y = y(:);

% Fitting a circle to the data - least squares style. 
%Here we start with
% (x-xo).^2 + (y-yo).^2 = R.^2
% Rewrite:
% x.^2 -2 x xo + xo^2 + y.^2 -2 y yo + yo.^2 = R.^2
% Put in matrix form:
% [-2x -2y 1 ] [xo yo -R^2+xo^2+yo^2]' = -(x.^2 + y.^2)
% Solve in least squares way...
A = [-2*x -2*y ones(length(x),1)];
x = A\-(x.^2+y.^2);
xo=x(1);
yo=x(2);
R = sqrt(  xo.^2 + yo.^2  - x(3));

end

function A = EllipseDirectFit(XY);
%
%  Direct ellipse fit, proposed in article
%    A. W. Fitzgibbon, M. Pilu, R. B. Fisher
%     "Direct Least Squares Fitting of Ellipses"
%     IEEE Trans. PAMI, Vol. 21, pages 476-480 (1999)
%
%  Our code is based on a numerically stable version
%  of this fit published by R. Halir and J. Flusser
%
%     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%
%     Output: A = [a b c d e f]' is the vector of algebraic 
%             parameters of the fitting ellipse:
%             ax^2 + bxy + cy^2 +dx + ey + f = 0
%             the vector A is normed, so that ||A||=1
%
%  This is a fast non-iterative ellipse fit.
%
%  It returns ellipses only, even if points are
%  better approximated by a hyperbola.
%  It is somewhat biased toward smaller ellipses.
%
centroid = mean(XY);   % the centroid of the data set

D1 = [(XY(:,1)-centroid(1)).^2, (XY(:,1)-centroid(1)).*(XY(:,2)-centroid(2)),...
      (XY(:,2)-centroid(2)).^2];
D2 = [XY(:,1)-centroid(1), XY(:,2)-centroid(2), ones(size(XY,1),1)];
S1 = D1'*D1;
S2 = D1'*D2;
S3 = D2'*D2;
T = -inv(S3)*S2';
M = S1 + S2*T;
M = [M(3,:)./2; -M(2,:); M(1,:)./2];
[evec,eval] = eig(M);
cond = 4*evec(1,:).*evec(3,:)-evec(2,:).^2;
A1 = evec(:,find(cond>0));
A = [A1; T*A1];
A4 = A(4)-2*A(1)*centroid(1)-A(2)*centroid(2);
A5 = A(5)-2*A(3)*centroid(2)-A(2)*centroid(1);
A6 = A(6)+A(1)*centroid(1)^2+A(3)*centroid(2)^2+...
     A(2)*centroid(1)*centroid(2)-A(4)*centroid(1)-A(5)*centroid(2);
A(4) = A4;  A(5) = A5;  A(6) = A6;
A = A/norm(A);

end  %  EllipseDirectFit

function [x,y,z] = SphereToCartesian(theta,phi,edgecoor_sph)
%given a point on the 2D surface map defined by phi, theta
%convert to cartesian coordinates
%use the edge coordinates in spherical coordinates to find closest point
%and retrieve the particle radius at that point, then convert to cartesian

%find distance (in theta/phi space) of each edge coordinate
dist_coors = sqrt((edgecoor_sph(:,1)-phi).^2 + (edgecoor_sph(:,2)-theta).^2);
%find closest coordinate
[closest_val, closest_index] = min(dist_coors);
r = edgecoor_sph(closest_index,3);

x = r * sin(phi) * cos(theta);
y = r * sin(phi) * sin(theta);
z = r * cos(phi);

end

function curvatures_withinShells = GetCurvatures_withinShells(MPStats_data,degranulation_data,shellsize,mc,gc)

curvatures_withinShells = struct('MC_Degran',[],'MC_Synapse',[],'GC_Degran',[],'GC_Synapse',[]);
% coords_cart = MPStats_data.edgecoor_cart;
% coords_sph = MPStats_data.edgecoor_sph(:,[1 2]);
% 2023-03-16: AHS identified that ^ these 2 lines need to be re-aligned by
% Daan's "Check Processing Status" function

    [coords_cart,~,coords_sph,~] = check_processing_status(MPStats_data);

    % Collect curvature data within differently-sized radial shells.
    for idx=1:size(degranulation_data.Spots_Degran,1)
        coords_spot = degranulation_data.Spots_Degran(idx,:);
        D_tospot_pol = sqrt((coords_spot(1)-coords_sph(:,1)).^2 + (coords_spot(2)-coords_sph(:,2)).^2);
        % I choose to replace coords_spot here, because we move from the best-estimate in polar space to the real coordinate in Cartesian space.
        coords_spot = coords_cart(D_tospot_pol==min(D_tospot_pol),:); 
        D_tospot_cart = sqrt((coords_spot(1)-coords_cart(:,1)).^2 + (coords_spot(2)-coords_cart(:,2)).^2 + (coords_spot(3)-coords_cart(:,3)).^2);
        % Filters for distance size to collect data within radial shells.
        for jdx=1:length(shellsize)
            idx_withinShell = find(D_tospot_cart < shellsize(jdx));
            mc_within = mc(idx_withinShell); gc_within = gc(idx_withinShell);
            curvatures_withinShells.MC_Degran{idx,jdx} = mc_within;
            curvatures_withinShells.GC_Degran{idx,jdx} = gc_within;
        end
    end

    % Collect curvature data within the entire synapse selection.
        roi_synapse = degranulation_data.Synapse_Coords_Sph;
        x = roi_synapse(:,1); y=roi_synapse(:,2);
        theta = coords_sph(:,1); phi = coords_sph(:,2);
        contact_map = inpolygon(theta,phi,x,y);     
        mc_within = mc(contact_map); gc_within = gc(contact_map);
        curvatures_withinShells.MC_Synapse{1,1} = mc_within;
        curvatures_withinShells.GC_Synapse{1,1} = gc_within;
    
end
