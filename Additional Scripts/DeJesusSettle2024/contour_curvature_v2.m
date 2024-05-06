%% Contour curvature analysis
% MDJ 20240126: procedural code to compile stacks from ImageJ
% CurvatureContour plugin

% v2: 2024-01-29: MH asked me to find a way to parse the data such that we
% can "unfold" the contours and "align" them to the mapped centroid

%% Initialization

path = uigetdir(pwd,'Select folder with .csvs of contour curvatures');
addpath(path);

%% Load matrix files into mother file for parsing & subsequent processing

% Note: readmatrix() generates a data array as:
%   index   x_coordinate  y_coordinate  x_contour   y_contour   curvature
%   1       2               3           4           5           6

count = 0;
contour_struct = struct('filename',{},'data',[],'unique_movie',{},'timepoint',[],'slice',[]);

file_list = dir(path);
for file = file_list'
    if contains(file.name,['.csv'])
        count = count+1;
        contour_struct(count).filename = file.name;
        contour_struct(count).data = readmatrix(file.name);
        contour_struct(count).unique_movie = extractBefore(...
            file.name,'_SegmentMask');
        contour_struct(count).timepoint = str2double(...
            extractBetween(file.name,'SegmentMask_','_'));
        contour_struct(count).slice = str2double(...
            extractBetween(file.name,'slice','_'));
    end
end

contour_struct = sortStruct(contour_struct,'slice');
contour_struct = sortStruct(contour_struct,'timepoint');
contour_struct = sortStruct(contour_struct,'unique_movie');

%% Parse the matrices and generate contour plots per frame

% Parsing is done by organizing the data first by cell, then by timepoint.
% (Slices should already be ordered)

% cell_idx is the master variable controlling which cell's timelapses to
% process at a time.

% Process by cell
list_cell = {contour_struct.unique_movie}';
unique_cell_ID = unique(list_cell);

msg = "Choose which cell to process:";
choice = menu(msg,unique_cell_ID);

for cell_idx=choice
    cell_ID = unique_cell_ID{cell_idx};
    [match_cell_ID] = ismember(list_cell,cell_ID);
    contour_bycell = contour_struct(match_cell_ID);
    
    % Within each cell, process by timepoint
    list_timepoint = [contour_bycell.timepoint];
    unique_time_ID = unique(list_timepoint);
    for jdx=1:length(unique_time_ID)
        time_ID = unique_time_ID(jdx);
        [match_time_ID] = ismember(list_timepoint,time_ID);
        contour_bycell_bytime = contour_bycell(match_time_ID);
        
        % Within each timepoint, process the stack:
        % Plot the image stack of contour curvatures
        figure, hold on
        for kdx=1:length(contour_bycell_bytime)
        contour_xyz = [...
            contour_bycell_bytime(kdx).data(:,4),...
            contour_bycell_bytime(kdx).data(:,5),...
            repmat(kdx,length(contour_bycell_bytime(kdx).data),1)];
        contour_curv = [contour_bycell_bytime(kdx).data(:,6)];
        scatter3(...
            contour_xyz(:,1),...
            contour_xyz(:,2),...
            contour_xyz(:,3),...
            5,contour_curv,'o','filled');
        end
        pbaspect([1 1 .1]),view(2),colormap(turquoisebrown),caxis([-0.05 0.05])
        title(contour_bycell(jdx).filename,'interpreter','none')
    end
    
end

clear contour_bycell contour_bycell_bytime contour_curv contour_xyz count idx jdx kdx...
    cell_idx choice match_cell_ID match_time_ID

%% Generate bounding box and compile data from in/out

roi_contours = struct('synapse_centroid',[],'synapse_roi',[],'synapse_xyz',[],...
    'synapse_curvature',[],'nonsynapse_centroid',[],'nonsynapse_roi',[],...
    'nonsynapse_xyz',[],'nonsynapse_curvature',[]);

% Find number of figures = timeframes of cell_ID
figlist = findobj('Type', 'figure');
n_frames = numel(figlist);

% Loop through each time frame and draw the polygons
for idx=1:n_frames
    
    figure(idx)
    disp(cell_ID)
    disp(idx)
    disp('Draw a region containing the synapse.')
       % Draw a polygon to define an ROI
        roi_selection = drawpolygon('FaceAlpha',0,'Color',[0 0 0],'StripeColor',[1 1 1],'LineWidth',1);
        % Subset the contour structure to the cell & appropriate frame
        [match_cell_ID] = ismember(list_cell,cell_ID);
        contour_bycell = contour_struct(match_cell_ID);
        [match_time_ID] = ismember(list_timepoint,idx);
        contour_bycell_bytime = contour_bycell(match_time_ID);
        % Compile all the slices into a big array; columns 4 & 5
        % contain the contour x and y, column 6 contains the curvature,
        % column 7 is appended with the slice number as a "z" value
        data_across_slices = [];
        for jdx=1:length(contour_bycell_bytime)
            data_across_slices = [data_across_slices;...
                contour_bycell_bytime(jdx).data,...
                repmat(...
                contour_bycell_bytime(jdx).slice,length(contour_bycell_bytime(idx).data),1)
                ];
        end
        % contour_xyz contains x-value, y-value, curvature, and z-slice
        contour_xyz = [data_across_slices(:,4),...
            data_across_slices(:,5),...
            data_across_slices(:,6),...
            data_across_slices(:,7)];
        within_roi = inpolygon(contour_xyz(:,1),contour_xyz(:,2),...
            roi_selection.Position(:,1),roi_selection.Position(:,2));
        data_within_roi = data_across_slices(within_roi,[4:7]);
        % Store info into the structure
        roi_contours(idx).synapse_roi = roi_selection.Position;
        roi_contours(idx).synapse_xyz = data_within_roi(:,[1 2 4]);
        roi_contours(idx).synapse_curvature = data_within_roi(:,3);
    disp('Mark the approximate cell front of the contact.')
        centroid = drawpoint('Color','r');
        roi_contours(idx).synapse_centroid = centroid.Position;
    disp('Draw a region outside the synapse as a control.')
       % Draw a polygon to define an ROI
        roi_selection = drawpolygon('FaceAlpha',0,'Color',[0 0 0],'StripeColor',[1 1 1],'LineWidth',1);
        % Subset the contour structure to the cell & appropriate frame
        [match_cell_ID] = ismember(list_cell,cell_ID);
        contour_bycell = contour_struct(match_cell_ID);
        [match_time_ID] = ismember(list_timepoint,idx);
        contour_bycell_bytime = contour_bycell(match_time_ID);
        % Compile all the slices into a big array; columns 4 & 5
        % contain the contour x and y, column 6 contains the curvature
%         data_across_slices_2 = [];
%         for kdx=1:length(contour_bycell_bytime)
%             data_across_slices_2 = [data_across_slices_2;contour_bycell_bytime(kdx).data];
%         end
        contour_xyz = [data_across_slices(:,4) data_across_slices(:,5)];
        within_roi = inpolygon(contour_xyz(:,1),contour_xyz(:,2),...
            roi_selection.Position(:,1),roi_selection.Position(:,2));
        data_within_roi = data_across_slices(within_roi,[4:7]);
        % Store info into the structure
        roi_contours(idx).nonsynapse_roi = roi_selection.Position;
        roi_contours(idx).nonsynapse_xyz = data_within_roi(:,[1 2 4]);
        roi_contours(idx).nonsynapse_curvature = data_within_roi(:,3);
    disp('Mark the centroid of the control region.')
        centroid = drawpoint('Color','g');
        roi_contours(idx).nonsynapse_centroid = centroid.Position;
end
   
%% Graph: Histogram comparing synapse vs. non-synapse curvature bins

figure
if exist('roi_contours','var')
    n_frames = length(roi_contours);
end

for idx=1:n_frames
    
    if n_frames > 5
        subplot(2,ceil(n_frames/2),idx),hold on
    elseif n_frames < 5
        subplot(1,n_frames,idx),hold on
    end

h_syn = histogram(roi_contours(idx).synapse_curvature,'BinWidth',.001,'Normalization','pdf','EdgeAlpha',0,'FaceColor',[.75 0 0]);
h_nonsyn = histogram(roi_contours(idx).nonsynapse_curvature,'BinWidth',.001,'Normalization','pdf','EdgeAlpha',0,'FaceColor',[.5 .5 .5]);
xlim([-0.05 0.05]), pbaspect([1 1 1])

end

%% Graph: manual norm2max histogram

figure
if exist('roi_contours','var')
    n_frames = length(roi_contours);
end

for idx=1:n_frames
    
    if n_frames > 5
        subplot(2,ceil(n_frames/2),idx),hold on
    elseif n_frames < 5
        subplot(1,n_frames,idx),hold on
    end
    
c_ns = roi_contours(idx).nonsynapse_curvature;
c_s = roi_contours(idx).synapse_curvature;
binW = 0.001;
[h_ns,edges] = histcounts(c_ns,'BinWidth',binW,'Normalization','count');
centers_ns = mean([edges(1:end-1);edges(2:end)]);
[h_s,edges] = histcounts(c_s,'BinWidth',binW,'Normalization','count');
centers_s = mean([edges(1:end-1);edges(2:end)]);
h_ns_scaled = h_ns/max(h_ns);
h_s_scaled = h_s/max(h_s);

% plot as bar plots with specified colors (or however you want to plot them)
area(centers_s,h_s_scaled,'EdgeAlpha',0,'FaceColor',[.75 0 0],'FaceAlpha',0.6)
area(centers_ns,h_ns_scaled,'EdgeAlpha',0,'FaceColor',[.5 .5 .5],'FaceAlpha',0.6)
xlim([-0.05 0.05]), pbaspect([1 1 1])

% histogram(roi_contours(idx).nonsynapse_curvature,'BinWidth',.001,'Normalization','probability','EdgeAlpha',0,'FaceColor',[.5 .5 .5])
% histogram(roi_contours(idx).synapse_curvature,'BinWidth',.001,'Normalization','probability','EdgeAlpha',0,'FaceColor',[.75 0 0])
% xlim([-0.05 0.05]), pbaspect([1 1 1])
end

%% Graph: Curvature vs. unfolded contour, per slice

% The data from the ImageJ plugin is already arranged according to
% next-neighbor connectivity, on account of them being sampled from a
% Fourier contour. So, to unfold and align, just pick the points in each
% slice nearest to the Centroid/Cell-Front-Point, then phase shift
% everybody such that the "center" point is at origin. Then, replace the
% notion of Euclidean pixel distance with topological/network/index
% distance from origin as new x-axis, and take curvatures as new y-axis.

figure, hold on
frames_to_analyze = 1:4;
for frame_idx = frames_to_analyze
    
    % Iterate over slices
    list_slices = unique(roi_contours(frame_idx).synapse_xyz(:,3));
    for idx = 1:length(list_slices)
        % Subset per slice
        [match_slice] = ismember(roi_contours(frame_idx).synapse_xyz(:,3),list_slices(idx));
        xyz_at_slice = roi_contours(frame_idx).synapse_xyz(match_slice,:);
        c_at_slice = roi_contours(frame_idx).synapse_curvature(match_slice,:);
        % Find point to anchor at origin
        centroid = roi_contours(frame_idx).synapse_centroid;
        d_to_centroid = sqrt((centroid(1) - xyz_at_slice(:,1)).^2 ...
            + (centroid(2) - xyz_at_slice(:,2)).^2);
        anchor_point = find(d_to_centroid == min(d_to_centroid));
        % Add topological reindexing with this point at origin
        idx_anchored = [1:length(d_to_centroid)] - anchor_point;
        subplot(2,length(frames_to_analyze),frame_idx), hold on
        scatter(idx_anchored,c_at_slice,5,'.');
        ylim([-.075 .075])
        xlim([-300 300])
        pbaspect([1 1 1])
        
        % Repeat for non-synapse regions
        [match_slice] = ismember(roi_contours(frame_idx).nonsynapse_xyz(:,3),list_slices(idx));
        xyz_at_slice = roi_contours(frame_idx).nonsynapse_xyz(match_slice,:);
        c_at_slice = roi_contours(frame_idx).nonsynapse_curvature(match_slice,:);
        % Find point to anchor at origin
        centroid = roi_contours(frame_idx).synapse_centroid;
        d_to_centroid = sqrt((centroid(1) - xyz_at_slice(:,1)).^2 ...
            + (centroid(2) - xyz_at_slice(:,2)).^2);
        anchor_point = find(d_to_centroid == min(d_to_centroid));
        % Add topological reindexing with this point at origin
        idx_anchored = [1:length(d_to_centroid)] - anchor_point;
        subplot(2,length(frames_to_analyze),frame_idx+length(frames_to_analyze)), hold on
        scatter(idx_anchored,c_at_slice,5,'.');
        ylim([-.075 .075])
        xlim([-300 300])
        pbaspect([1 1 1])
    end
    
end

%% Graph: Curvature vs. unfolded contour MEAN LINES 

plot_with_errors = 0;

figure, hold on
frames_to_analyze = 1:length(roi_contours);
for frame_idx = frames_to_analyze
    
    curvatures_by_slice = struct('synapse_topoidx',[],'synapse_curvature',[],...
        'nonsynapse_topoidx',[],'nonsynapse_curvature',[]);    
    % Iterate over slices
    list_slices = unique(roi_contours(frame_idx).synapse_xyz(:,3));
    for idx = 1:length(list_slices)
        % Subset per slice
        [match_slice] = ismember(roi_contours(frame_idx).synapse_xyz(:,3),list_slices(idx));
        xyz_at_slice = roi_contours(frame_idx).synapse_xyz(match_slice,:);
        c_at_slice = roi_contours(frame_idx).synapse_curvature(match_slice,:);
        % Find point to anchor at origin
        centroid = roi_contours(frame_idx).synapse_centroid;
        d_to_centroid = sqrt((centroid(1) - xyz_at_slice(:,1)).^2 ...
            + (centroid(2) - xyz_at_slice(:,2)).^2);
        anchor_point = find(d_to_centroid == min(d_to_centroid),1,'first');
        % Add topological reindexing with this point at origin
        idx_anchored = [1:length(d_to_centroid)] - anchor_point;
        curvatures_by_slice(idx).synapse_topoidx = idx_anchored;
        curvatures_by_slice(idx).synapse_curvature = c_at_slice;
        
        % Repeat for non-synapse regions
        [match_slice] = ismember(roi_contours(frame_idx).nonsynapse_xyz(:,3),list_slices(idx));
        xyz_at_slice = roi_contours(frame_idx).nonsynapse_xyz(match_slice,:);
        c_at_slice = roi_contours(frame_idx).nonsynapse_curvature(match_slice,:);
        % Find point to anchor at origin
        centroid = roi_contours(frame_idx).nonsynapse_centroid;
        d_to_centroid = sqrt((centroid(1) - xyz_at_slice(:,1)).^2 ...
            + (centroid(2) - xyz_at_slice(:,2)).^2);
        anchor_point = find(d_to_centroid == min(d_to_centroid),1,'first');
        % Add topological reindexing with this point at origin
        idx_anchored = [1:length(d_to_centroid)] - anchor_point;
        curvatures_by_slice(idx).nonsynapse_topoidx = idx_anchored;
        curvatures_by_slice(idx).nonsynapse_curvature = c_at_slice;
    end
    
    % Generating meanlines plots:
    
    idx_synapse_aligned = nan(length(curvatures_by_slice),20001);
    c_synapse_aligned = nan(length(curvatures_by_slice),20001);
    for idx=1:length(curvatures_by_slice)
        % These odd conditionals are necessary to exclude extreme slices where the
        % cell is completely cut off from the ROI or doesn't exist at that
        % height anymore.
        if size(curvatures_by_slice(idx).synapse_curvature)>0
            i_syn = curvatures_by_slice(idx).synapse_topoidx;
            c = curvatures_by_slice(idx).synapse_curvature;
            idx_synapse_aligned(idx,10001+i_syn(1):10001+i_syn(end)) = i_syn;
            c_synapse_aligned(idx,10001+i_syn(1):10001+i_syn(end)) = c;
        end
    end
    idx_nonsynapse_aligned = nan(length(curvatures_by_slice),20001);
    c_nonsynapse_aligned = nan(length(curvatures_by_slice),20001);
    for idx=1:length(curvatures_by_slice)
        if size(curvatures_by_slice(idx).nonsynapse_curvature)>0
            i_nonsyn = curvatures_by_slice(idx).nonsynapse_topoidx;
            c = curvatures_by_slice(idx).nonsynapse_curvature;
            idx_nonsynapse_aligned(idx,10001+i_nonsyn(1):10001+i_nonsyn(end)) = i_nonsyn;
            c_nonsynapse_aligned(idx,10001+i_nonsyn(1):10001+i_nonsyn(end)) = c;
        end
    end
    
    % For plotting: squish (average down rows) indices & curvatures
    x_nonsyn = mean(idx_nonsynapse_aligned,1,"omitnan");
    y_nonsyn = mean(c_nonsynapse_aligned,1,"omitnan");
    sd_nonsyn = std(c_nonsynapse_aligned,1,"omitnan");    
    x_syn = mean(idx_synapse_aligned,1,"omitnan");
    y_syn = mean(c_synapse_aligned,1,"omitnan");
    sd_syn = std(c_synapse_aligned,1,"omitnan");    

    if length(roi_contours) > 5
        subplot(2,length(frames_to_analyze)/2,frame_idx), hold on
    else
        subplot(1,length(frames_to_analyze),frame_idx),hold on
    end
    if plot_with_errors == 1
        leftbound = find(~isnan(x_nonsyn),1,'first');
        rightbound = find(~isnan(x_nonsyn),1,'last');
        fill([x_nonsyn(leftbound:rightbound),fliplr(x_nonsyn(leftbound:rightbound))],...
            [y_nonsyn(leftbound:rightbound)-sd_nonsyn(leftbound:rightbound),flipud(y_nonsyn(leftbound:rightbound)+sd_nonsyn(leftbound:rightbound))],...
            [.75 .75 .75],'linestyle','none','FaceAlpha',0.5)
    end
    plot(x_nonsyn,y_nonsyn,'Color',[.5 .5 .5],'LineWidth',2.5)
    if plot_with_errors == 1
        leftbound = find(~isnan(x_syn),1,'first');
        rightbound = find(~isnan(x_syn),1,'last');
        fill([x_syn(leftbound:rightbound),fliplr(x_syn(leftbound:rightbound))],...
            [y_syn(leftbound:rightbound)-sd_syn(leftbound:rightbound),flipud(y_syn(leftbound:rightbound)+sd_syn(leftbound:rightbound))],...
            [.75 .5 .5],'linestyle','none','FaceAlpha',0.5)
    end
    plot(x_syn,y_syn,'Color',[.75 0 0],'LineWidth',2.5)
    
    xlim([-500 500])
    ylim([-0.025 0.025])
    pbaspect([1 1 1])
     
end
