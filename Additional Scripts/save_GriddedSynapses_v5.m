%% saveGriddedSynapses.m (V5)

%STARTED 02/01/2021
%This script is written to take a data_structure of MPRenders and ROIs (as
%done at the start of stain_correlation_plots) and create a new
%data structure containing interpolated NxN grids encoding curvature and
%stain data (if applicable). 

% 2021-02-12...15: MDJ edited to find the maximum distance from centroid,
% through a curvature point of interest, up to the cell edge. This is so we
% can re-define our radial distance = 1.0 as the centroid-to-edge distance
% at that angle.

% 2021-03-09: MDJ edited with remarks for non-stained datasets.
% 2021-07-15: MDJ edited to combine with Gaussian curvatures.
% 2021-08-05: Added actin intensity gridding and distance_data structure.

% 2021-08-05: Unified MeanC, GaussC, and Cellstain in the distance_data structure. 

%% Initialize or Load Actin/Curvature Data Structure
% Check for existing data structure. Load it into memory if it exists,
% create a new structure if it doesn't. 

clear('struct_file','struct_path')
stain_check = 2;
    while ~ismember(stain_check,[0 1])
        stain_check = input("Include cell stain data? 1=Yes, 0=No ");
    end

name_format = input("How do you want to name the data output structure? ",'s');
file_struct_name = strcat(name_format,'.mat');

current_dir = pwd;
[struct_file,struct_path] = uigetfile('*.mat',"Select existing data file, hit cancel if none exists");

if struct_file
    cd(struct_path)
    load(struct_file)
else
    data_struct = struct('Label',{},'MPRender_label',{},'ROI_label',{},'edgecoor',{},'edgecoor_contact',{}, ...
        'edgecoor_sph',{},'edgecoor_sph_contact',{},'mc',{},'mc_contact',{},'gc',{}, ...
        'gc_contact',{},'stain_int',{}, ...
        'stain_int_contact',{});
    if stain_check==0
        data_struct = rmfield(data_struct,'stain_int');
        data_struct = rmfield(data_struct,'stain_int_contact');
    end
    save(file_struct_name,'data_struct')
end

%% Select Folder of interest
% Define where the MPRender and ROI files are located (must be in same
% folder)
file_path = uigetdir(pwd,'Select folder containing MPRenders and ROI data');

%% Create or extend data structure (only do first time, or if adding more files)
method = 'ROI';

current_dir = pwd;
cd(file_path);
files = dir(file_path);

%Check current size of the data structure
current_data_size = length(data_struct);
k=1+current_data_size;

invalid_files = struct();
invalid_file_count = 0;

for file = files'
    if contains(file.name,'MPRender')
        fprintf('Loading %s \n',file.name)
        
        load(file.name)
        
        %check size of MPStats
        mp_size = length(MPStats);
        fprintf(['Number of Frames: ', num2str(mp_size),' \n'])
        
        %check if MPStats is already in dataframe
        current_labels = [];
        for n = 1:length(data_struct)
            current_labels{n} = getfield(data_struct(n),'Label');
        end
        match_check = strcmp(current_labels,MPStats(1).FileName);
        if any(match_check)
            [~,indx] = max(match_check);
            fprintf(['MPStats already in data_struct, starting at ' ...
                'index %d \n \n \n '], indx)
            continue
        end
        
        %Check if ROI file exists
        ROI_filename = [file.folder,'/', strrep(file.name,'MPRender','ROI')];
        if isfile(ROI_filename)
            load(ROI_filename)
        else
            disp('Error: No ROI file detected for...')
            roi_position = struct();
            disp(file.name)
        end               
        
        for j = 1:length(MPStats)
            %check if data is usable (various filters can be added here)
             if stain_check==1
                 valid = 1;
                  if ~isfield(MPStats(j),'isincontact')
                      valid = 0;
                      reason = 'isincontact field does not exist.';
                  elseif isempty(MPStats(j).isincontact)
                      valid = 0;
                      reason = 'isincontact is empty.';
                  end
             end
             
      % Comment out the lines between "valid = 1" and "end" above if no stain data exists.
            
            %If Invalid, store metadata 
            if valid == 0
                invalid_file_count = invalid_file_count+1;
                invalid_files(invalid_file_count).Label = MPStats(j).FileName;
                invalid_files(invalid_file_count).Reason_for_Error = reason;
                continue
            end            

            fprintf('Extracting data for %s \n',MPStats(j).FileName);
            fprintf('data_struct Index: %d \n \n \n', k) 
            new_stats = extract_stain_curv_data(MPStats,j,method,roi_position);
            data_struct(k) = new_stats;
            data_struct(k).Label = MPStats(j).FileName;
            
            % Added stuff for making ROI matching easier
            data_struct(k).MPRender_label = file.name;
            data_struct(k).ROI_label = ROI_filename;
            
            k=k+1;
            fprintf([num2str(k),'\n'])
        end
    end
end

cd(current_dir)
fprintf('Saving updated data_struct \n \n \n')
save(file_struct_name,'data_struct')

%% Interpolate and store gridded data: Mean curvatures, Gaussian curvatures, Cell Stain Intensity

out_path = uigetdir('pwd','Select output folder for grids');
cd(out_path)

for k = 1:length(data_struct)
    mc = data_struct(k).mc_contact;
    gc = data_struct(k).gc_contact;
    edgecoor = data_struct(k).edgecoor_sph_contact;
    if stain_check==1
        cell_stain = data_struct(k).stain_int_contact;
    end
    
    [~,~,mc_i] = ... 
        interpolateCurvature(edgecoor(:,1), edgecoor(:,2), ...
        mc,50);
    mc_i(isnan(mc_i)) = 0;
    save_name = strrep(data_struct(k).Label,'.tif','_50x50_grid_MeanC.mat');
    fprintf('Saving %s \n',save_name)
    save(save_name,'mc_i');
    
    [~,~,gc_i] = ...
        interpolateCurvature(edgecoor(:,1), edgecoor(:,2),gc,50);
    gc_i(isnan(gc_i)) = 0;
    save_name = strrep(data_struct(k).Label,'.tif','_50x50_grid_GaussC.mat');
    fprintf('Saving %s \n',save_name)
    save(save_name,'gc_i'); 
    
    % Comment out this segment if you are avoiding Stain Intensity grid
    % data.
    [~,~,int_i] = ... 
        interpolateCurvature(edgecoor(:,1), edgecoor(:,2), ...
        double(cell_stain),50);
    int_i(isnan(int_i)) = 0;
    save_name = strrep(data_struct(k).Label,'.tif','_50x50_grid_Cellstain.mat');
    fprintf('Saving %s \n',save_name)
    save(save_name,'int_i');
    
end 

fprintf('Done saving gridded interpolations.')
    
%% Calculate curvature/distance data

out_path = uigetdir('pwd','Select output folder for radial distance data.');
cd(out_path)

for k = 1:length(data_struct)
    edgecoor = data_struct(k).edgecoor_sph_contact;
    mc = data_struct(k).mc_contact;
    gc = data_struct(k).gc_contact;
    if stain_check==1
        cell_stain = data_struct(k).stain_int_contact;
    end    
    centroid =[mean(edgecoor(:,1)),mean(edgecoor(:,2))];
    distances = sqrt((edgecoor(:,1)'-centroid(1)).^2 + (edgecoor(:,2)'-centroid(2)).^2);
    
    % MDJ: Inserted here my function to obtain distances to the cell edge.
        % 1. Find the right ROI frame: Loop through the ROI labels until finding the matching one with data_struct(k).Label      
            load(data_struct(k).ROI_label);
            
            for idx=1:length(roi_position)
                if strcmp(roi_position(idx).label,data_struct(k).Label)==1
                    roi_idx=idx;
                end
            end
            
        % 2. Execute computations per curvature value
            roi_coords = roi_position(roi_idx).position;
            roi_coords(end+1,:) = roi_coords(1,:);
            roi_shape = polyshape(roi_coords);
            intersections = [];
            for idx=1:length(edgecoor)
                endoflinex=centroid(1); endofliney=centroid(2); p=1;
                while isinterior(roi_shape,endoflinex,endofliney)
                    endoflinex = centroid(1) + (edgecoor(idx,1)-centroid(1))*10^(p);
                    endofliney = centroid(2) + (edgecoor(idx,2)-centroid(2))*10^(p);
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
            
        if stain_check==1
    distance_data = struct('DistanceFromCentroid',distances,...
        'MeanCurvature',mc,'GaussianCurvature',gc,'DistanceFromEdge',distances_tocelledge','NormalizedRadialD',d_norm','StainInt',cell_stain');
        elseif stain_check==0
    distance_data = struct('DistanceFromCentroid',distances,...
        'MeanCurvature',mc,'GaussianCurvature',gc,'DistanceFromEdge',distances_tocelledge','NormalizedRadialD',d_norm');
        end
        
    % Removed StainInt',stain' from distance_data initialization
    
      save_name = strrep(data_struct(k).Label,'.tif','_distance_data.mat');
      fprintf('Saving %s \n',save_name)
      save(save_name,'distance_data');

end

fprintf('Done saving distance_data structures.')
    
%% Functions for setting plot options

% Set 3D plots options

function[ph] = Make_3D_plot(TRI,vert,color)

    ph = patch('faces',TRI,'vertices',vert,'facevertexcdata',mean(color(TRI),2),'facecolor','flat','edgecolor','none');
    view(3)
    axis equal
    zoom(1.4)
    axis off
    set(gca,'CameraViewAngleMode','Manual','Clipping','off')
    view(0,0)

end

% Set 2D plot options (by surface)
function[ph] = Make_2D_plot(TRI,vert,color)

    ph = patch('faces',TRI,'vertices',vert,'facevertexcdata',mean(color(TRI),2),'facecolor','flat','edgecolor','none');
    view(3)
    axis off
    axis equal
    view(0,90);
    set(gca,'Clipping','off')

end

% Set 2D plot options (by point)
function[ph] = Make_2D_dotplot(edges,color)
    scatter3(edges(:,1),edges(:,2),edges(:,3),10,color,'filled')
    
end

function[struct_line] = extract_stain_curv_data(MPStats,iMP,method,roi_position)
    struct_line = struct();
    struct_line.Label={};   
    struct_line.MPRender_label = {};
    struct_line.ROI_label = {};
    
    [edgecoor,centroidcoor,edgecoor_sph,TRI_nosph] = check_processing_status(MPStats(iMP));
    if strcmp(method,'isincontact')
        contact_map = MPStats(iMP).isincontact;
    elseif strcmp(method,'ROI')
        theta = edgecoor_sph(:,1); phi = edgecoor_sph(:,2);
        x = roi_position(iMP).position(:,1); y = roi_position(iMP).position(:,2); 
        contact_map = inpolygon(theta,phi,x,y); 
    else
        fprintf('Error, method for selecting contact area not recognized')
        return
    end
    
    ThetaPhiR = edgecoor_sph;
    R = ThetaPhiR(:,3);
    R_contact = R(contact_map);
    
    if isfield(MPStats,'stain_int')
        if strcmp(MPStats(iMP).Gridding,'equi')
            stainint    = MPStats(iMP).stain_int;
            stainintint = MPStats(iMP).stain_int_integrated;
        elseif isfield(MPStats,'stain_int_irr')
                    stainint    = MPStats(iMP).stain_int_irr;
                    stainintint = MPStats(iMP).stain_int_integrated_irr;
        else
            stainint    = MPStats(iMP).stain_int;
            stainintint = MPStats(iMP).stain_int_integrated;
            [~,stainint,~,stainintint] = ...
                Determine_Coverage(imgaussfilt3(MPStats(iMP).IMstain,1), MPStats(iMP).Weighted_Centroid'/MPStats(iMP).PixelSizeXY,...
                MPStats(iMP).EquivDiameter,MPStats(iMP).edgecoor_sph_wcatorigin(:,1),MPStats(iMP).edgecoor_sph_wcatorigin(:,2),[MPStats(iMP).PixelSizeXY MPStats(iMP).PixelSizeZ],...
                'R',MPStats(iMP).edgecoor_sph_wcatorigin(:,3)/MPStats(iMP).PixelSizeXY,'Gaussiansignal',MPStats(iMP).IMstain_isgaussian,'radialinterval' ,...
                MPStats(iMP).IMstain_radialinterval);
        end
        
        
    end
    
    [ThetaPhiRsmooth] = Smooth_spherical_mesh(edgecoor_sph,sqrt(1^2/pi),'arc',1,'mean');
    [X,Y,Z] = sph2cart(ThetaPhiRsmooth(:,1),ThetaPhiRsmooth(:,2),ThetaPhiRsmooth(:,3));    

    FV.faces    = MPStats(iMP).TRI_Connectivity_SPHbound;
    FV.vertices = [X,Y,Z];
    pc = GetCurvatures(FV,0);
    gc = pc(1,:).*pc(2,:);
    mc = mean(pc);
    
    struct_line.edgecoor = edgecoor;    
    struct_line.edgecoor_sph = edgecoor_sph;
    struct_line.edgecoor_contact = edgecoor(contact_map,:);
    struct_line.edgecoor_sph_contact = edgecoor_sph(contact_map,:);
    struct_line.mc = mc;
    struct_line.gc = gc;
    struct_line.mc_contact = mc(contact_map);
    struct_line.gc_contact = gc(contact_map);
    % Comment out the bottom lines if no stain data.
    if isfield(MPStats,'stain_int')
        struct_line.stain_int = stainint;
        struct_line.stain_int_contact = struct_line.stain_int(contact_map);  
    end

end

function [t_i,p_i,mc_i] = interpolateCurvature(theta,phi,mc,N)
%Given theta-phi coordinates of a sphere with mean curavture values
%associated with each point, create a regular NxN grid, using the MATLAB
%griddata function to calculate a rolling average for each grid square. The
%grid lines are determined by max/min theta-phi coordinates equally spaced.

t_grid = linspace(min(theta),max(theta),N);
p_grid = linspace(min(phi),max(phi),N);

[t_i,p_i] = meshgrid(t_grid,p_grid);

mc_i = griddata(theta,phi,mc,t_i,p_i);

end
