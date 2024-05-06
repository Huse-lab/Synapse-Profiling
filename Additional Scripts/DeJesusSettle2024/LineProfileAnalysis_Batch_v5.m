%Batch Processing for Line Profile Analysis of deformed microparticles
% Alex Settle, 10/18/22, v1 
% Put matching MPRender and ROI pools in the folder of interest and run
% step 1, followed by step 2 to initialize. Step 3 will take the longest,
% and create a new directory "Line_Profiles" within the existing directory
% and place the output CSVs and Matlab structures (LPA_Stats), which
% contain the traces. No plotting or individual line_trace exports yet. 
% Requirements: MATLAB 2021b and the circle_fit function (see within the
% folder)

% Edited graphing portion by MDJ 2022-10-18 for IMP Retreat talk.
% Edited graphing portion by MDJ 2022-10-26 for generating and saving
% V3 - major restructure to do 3D slicing by AHS 2023-01-31
% V4 - bug fix to get theta bounds at each slice to help with plotting
% V5 - bug fix to calculate Hertz force
% figures individually

% Note as of 2023-02-01 (MDJ): This code does not calculate the Hertz force
% estimate in the preferred way! We should refer to the old LPA_v2.m data for those
% violins.

%% Initialization: compile MPRender and ROI data into a single structure for analysis

LPA_Stats3D = struct('FileName',{},'edgecoor_cart_aligned',[],'edgecoor_sph_aligned',[],'stain_int',[],'roi_position',[],'ROI_Z_bounds',[],'ROI_Theta_Bounds',[],'Num_ZSlices',[],...
        'LineTraceTables',[],'CenterDepth_um',[],'Center_Slice_Deviation_ROI',[],'Center_Slice_Deviation_Total',[], ...
        'Volume_Deviation_ROI',[],'Volume_Deviation_Total',[],'Est_Force_Hertz_pN',[],...
        'Est_BulkStress_Pa',[],'Est_BulkForce_pN',[],'ParticleRadius',[],'InitialParticleVolume',[],...
        'dVStrain',[]);

path = uigetdir(pwd,'Select folder containing MPRender and ROI data for indentation analysis.');
figlist = dir(fullfile(path, '*MPRender_pool*'));

% Loop to fetch data from each data structure and stuff them into the
% mother structure.
for idx=1:length(figlist)
    % Load MPRender and ROI files.
    MPname = fullfile(path,figlist(idx).name);
    load(MPname);
    disp(figlist(idx).name);
    ROIname = strrep(MPname,'MPRender_pool','ROI_pool');
    load(ROIname);
    % Stuff data row by row into the mother structure for calculation.
    for jdx=1:length(MPStats)
        counter = length(LPA_Stats3D)+1;
        LPA_Stats3D(counter).FileName = MPStats(jdx).FileName;
        LPA_Stats3D(counter).edgecoor_cart_aligned = MPStats(jdx).edgecoor_cart_aligned;
        LPA_Stats3D(counter).edgecoor_sph_aligned = MPStats(jdx).edgecoor_sph_aligned;
        LPA_Stats3D(counter).stain_int = MPStats(jdx).stain_int;
        LPA_Stats3D(counter).roi_position = roi_position(jdx).position;
    end
end

%% Calculate volume profiles (LPA_Stats3D)
% This is the last code block updated by Alex on 1/31/23
% Beyond this is the same as v2
num_slices = 20;

for n = 1:length(LPA_Stats3D)
    try
        volumeStats = calculate_volumeStats(LPA_Stats3D(n),num_slices);
        LPA_Stats3D(n) = volumeStats;
    catch
        fprintf(['Failure on ' LPA_Stats3D(n).FileName ', leaving empty \n'])
    end
end

LPA_Stats3D = assign_celltype(LPA_Stats3D);
LPA_Stats3D = sortStruct(LPA_Stats3D,'Celltype');

%% Visualize indentation traces one by one

% Assign idx
k = 1:length(LPA_Stats3D);
idx = 10;

    figure
    plot_blobtrace(LPA_Stats3D(k(idx)));
    figure
    plot_lineprofile(LPA_Stats3D(k(idx)));

%% Plot and save all traces

k = 1:length(LPA_Stats3D);

% Plot all BlobTraces first

    for idx=1:length(k)
        figure
        plot_blobtrace(LPA_Stats3D(k(idx)));
        saveas(gcf,strcat('BlobTrace_Frame',num2str(k(idx)),'.tif'))
    end
    close all

% Plot all LineProfiles, then standardize their xlim & ylim values

    xbounds = []; ybounds = [];
    for idx=1:length(k)
        figure
        plot_lineprofile(LPA_Stats3D(k(idx)));
        xbounds(idx,:) = xlim; ybounds(idx,:) = ylim;
    end

% Standardize xlim & ylim across all open figures

    % Find max ranges for x & y, and apply them to everybody
    xbounds_forAll = [min(xbounds(:)) max(xbounds(:))];
    ybounds_forAll = [min(ybounds(:)) max(ybounds(:))];
    % Save and close loop
    for idx=1:length(LPA_Stats3D)
        figure(idx)
        xlim(xbounds_forAll); ylim(ybounds_forAll);
        saveas(gcf,strcat('LineProfile_Frame',num2str(k(idx)),'.tif'));
    end
    close all

%% Functions

% 3D Stats function
function volumeStats = calculate_volumeStats(lpa_input,num_slices)

    C = lpa_input.edgecoor_cart_aligned;
    C_sphere = lpa_input.edgecoor_sph_aligned;
    stain = lpa_input.stain_int;
    R = mean(C_sphere(:,3));
    ROImask = inpolygon(C_sphere(:,1),C_sphere(:,2),...
        lpa_input.roi_position(:,1),...
        lpa_input.roi_position(:,2));
    %plot_blobtrace(LPA_Stats)
    z_range = [min(C(ROImask,3)),max(C(ROImask,3))];
    z_sections = linspace(z_range(1),z_range(2),num_slices+1);
    z_height = range(z_range);
    center_range = mean(z_range)+[-z_height*.15 z_height*.15];

    sectionStats = [];
    
    for k = 1:length(z_sections)-1
        slice_Stats = calculate_Sliceprofile(lpa_input,[z_sections(k) z_sections(k+1)],5);
        sectionStats = [sectionStats slice_Stats];
    end

    centerslice_Stats = calculate_Sliceprofile(lpa_input,center_range,50);

    volumeStats = lpa_input;
    volumeStats.ROI_Z_bounds = z_range;
    volumeStats.Num_ZSlices = num_slices;
    volumeStats.LineTraceTables = {sectionStats(:).LineTraceTable};

    volumeStats.ROI_Theta_Bounds = {sectionStats(:).ROI_Theta_Bounds};

    volumeStats.ParticleRadius = max([sectionStats(:).IdealRadius]);
    volumeStats.InitialParticleVolume = (4/3)*pi*volumeStats.ParticleRadius^3;

    volDevROI = trapz(z_sections(1:end-1),...
        [sectionStats(:).Integrated_Deviation_ROI]);
    volDevTotal = trapz(z_sections(1:end-1),...
        [sectionStats(:).Integrated_Deviation_Total]);
    volumeStats.Volume_Deviation_ROI = volDevROI;
    volumeStats.Volume_Deviation_Total = volDevTotal;

    

    volumeStats.Est_Force_Hertz_pN = centerslice_Stats.Estimated_Force_pN;
    volumeStats.dVStrain = abs(volumeStats.Volume_Deviation_ROI / ...
        volumeStats.InitialParticleVolume);
    
    Eparticle = 300; %N/m^2 (Pa) 
    volumeStats.Est_BulkStress_Pa = Eparticle * volumeStats.dVStrain;
    est_surfaceArea = pi*((z_range(2)-z_range(1))/2)^2; %µm^2
    est_totalForce_N = est_surfaceArea*volumeStats.Est_BulkStress_Pa*10^-12;
    volumeStats.Est_BulkForce_pN = est_totalForce_N*10^12;
end

% Adjusted form of calculate_lineprofile to be compatible with 3D analysis
function slice_Stats = calculate_Sliceprofile(lpa_input,z_band,smoothing_factor)

        C = lpa_input(1).edgecoor_cart_aligned;
        C_sphere = lpa_input(1).edgecoor_sph_aligned;
        stain = lpa_input(1).stain_int;
        ROImask = inpolygon(C_sphere(:,1),C_sphere(:,2),...
            lpa_input(1).roi_position(:,1),lpa_input(1).roi_position(:,2));
        
        % Get center slice based on ROI center in spherical and cartesian
        % Project slice onto 2D-space and get ROI mask in 2D
        band_mask = C(:,3) > z_band(1) & C(:,3) < z_band(2);
        Csphere_filt = C_sphere(band_mask,:);
        stain_filt = stain(band_mask);
        Cfilt = C(band_mask,:);
        roifilt = ROImask(band_mask);
        
        
        % Fit Circle to non-ROI points
        [xo,yo,R0] = circle_fit(Cfilt(~roifilt,1),Cfilt(~roifilt,2));
        zo = mean(Cfilt(:,3));

        % Convert 2D cartesian data to spherical with new center (subtract center)
        Cshift = Cfilt-[xo,yo,zo];
        [theta,rho] = cart2pol(Cshift(:,1),Cshift(:,2));
        Cshift_pol = [theta,rho];
        % Calculate rolling average to smooth line profile
        theta_space = linspace(-pi,pi,500);
        [sortedC,sort_idx] = sortrows(Cshift_pol);
        sortedROI = roifilt(sort_idx);
        rollingav = movmean(sortedC(:,2),smoothing_factor);
        rollingav_function = interp1(sortedC(:,1),rollingav,theta_space);
        ROI_function = interp1(sortedC(:,1),double(sortedROI),theta_space);
        rollingav_function(isnan(rollingav_function)) = R0;
        ROI_function(isnan(ROI_function))= 0;
        
        %Also calculate the flat function at R0;
        flat_function = ones(length(theta_space),1)*R0;

        ROIbound1 = find(ROI_function,1,'first');
        ROIbound2 = find(ROI_function,1,'last');

        %Calculate integral of deviation from ideal
    % This is the original line of code, which should be applied for T vs.
    % T synapses... It integrates positive and negative deviations
    % together.
%         deviation_function = abs(flat_function-rollingav_function');
    % This is the alternative line of code for emphasizing the
    % directionality differences in deviation from T cells vs. phagocytes.
        deviation_function = flat_function-rollingav_function';
        
        devIntegral = trapz(theta_space*R0,deviation_function);
        devIntegralROI = trapz(R0*theta_space(ROIbound1:ROIbound2),...
            deviation_function(ROIbound1:ROIbound2));

        %Estimate Force cell assumed to be a rigid half-plane
        Eparticle = 300; %N/m^2 (Pa)
        particle_radius = mean(flat_function);
        min_radius = min(rollingav_function);
        d = (particle_radius - min_radius)*10^-6; %m
        v = 0.457; % get citation. poisson ratio of pAAm
        cell_radius = inf; % µm (assumption kept constant across all cells)
        Eprime = Eparticle / (1-v^2);
        Reff = (1/ (1/particle_radius + 1/cell_radius))*(10^-6); %m
        F = (4/3)*Eprime*(Reff^(1/2))*(d^(3/2)); % Newton
        FpN = F * 10^12; %pico-Newton

        slice_Stats(1).ROI_Theta_Bounds = [theta_space(ROIbound1),theta_space(ROIbound2)];
        slice_Stats(1).Z_Band = [z_band(1),z_band(2)];
        slice_Stats(1).LineTraceTable = table(theta_space',theta_space'*2*R0,...
            rollingav_function',flat_function,'VariableNames',{'Theta',...
            'Perimeter_Trace','Measured_Radius','Idealized_Radius'});
        slice_Stats(1).MaximumDepth_um = d*10^6;
        slice_Stats(1).IdealRadius = R0;
        slice_Stats(1).Integrated_Deviation_Total  = devIntegral;
        slice_Stats(1).Integrated_Deviation_ROI  = devIntegralROI;
        slice_Stats(1).Estimated_Force_pN = FpN;
        
%     save([out_path '/' outFile],'LPA_Stats')
%     outcsv = strrep(outFile,'.mat','.csv');
%     outCSVTable = table( (1:length(LPA_Stats))',[LPA_Stats.MaximumDepth_um]',...
%         [LPA_Stats.Integrated_Deviation_ROI]',...
%         [LPA_Stats.Integrated_Deviation_Total]',...
%         'VariableNames',{'Index','Crater Depth (µm)','Integrated Deviation on ROI (µm^2)', ...
%         'Integrated Deviation Total'});
%     writetable(outCSVTable,outcsv)

end

function plot_blobtrace(stats)

    hold on
    C = stats(1).edgecoor_cart_aligned;
    C_sphere = stats(1).edgecoor_sph_aligned;
    stain = stats(1).stain_int;
    R = mean(C_sphere(:,3));
    ROImask = inpolygon(C_sphere(:,1),C_sphere(:,2),...
        stats(1).roi_position(:,1),stats(1).roi_position(:,2));

    Cfilt = C(abs(C_sphere(:,2)) < 0.2,:);
    roifilt = ROImask(abs(C_sphere(:,2)) < 0.2);

    % Add "expected radius" for the ROI region
    [xo,yo,R] = circle_fit(Cfilt(~roifilt,1),Cfilt(~roifilt,2));
    fitted_theta = linspace(stats(1).ROI_Z_bounds(1),stats(1).ROI_Z_bounds(2),70);
    fitted_x = R*cos(fitted_theta) + xo;
    fitted_y = R*sin(fitted_theta) + yo;
    scatter(fitted_x,fitted_y,5,'k','o','filled');
 
%Color assignment
kulay = [0 0 1];
    scatter(Cfilt(roifilt,1),Cfilt(roifilt,2),5,kulay,'o','filled');
    kulay = [.75 .75 .75];
    scatter(Cfilt(~roifilt,1),Cfilt(~roifilt,2),5,kulay,'o','filled');    
    axis off
    pbaspect([1 1 1]), daspect([1 1 1])

end

function plot_lineprofile(stats)

    hold on
    slice_to_plot = ceil(size(stats(1).LineTraceTables,2)/2);
    roi_plot_range = stats(1).LineTraceTables{1,slice_to_plot}.Theta > stats(1).ROI_Theta_Bounds{1,slice_to_plot}(1)...
        & stats(1).LineTraceTables{1,slice_to_plot}.Theta < stats(1).ROI_Theta_Bounds{1,slice_to_plot}(2);
    plot(stats(1).LineTraceTables{1,slice_to_plot}.Perimeter_Trace,...
        stats(1).LineTraceTables{1,slice_to_plot}.Idealized_Radius,'k:','LineWidth',1)
    plot(stats(1).LineTraceTables{1,slice_to_plot}.Perimeter_Trace,...
        stats(1).LineTraceTables{1,slice_to_plot}.Measured_Radius,'Color',[.75 .75 .75],'LineWidth',3)
    plot(stats(1).LineTraceTables{1,slice_to_plot}.Perimeter_Trace(roi_plot_range),...
        stats(1).LineTraceTables{1,slice_to_plot}.Measured_Radius(roi_plot_range),'b-','LineWidth',3);
    % Set ylim based on idealized radius.
    yband = ceil(abs(min(stats(1).LineTraceTables{1,slice_to_plot}.Measured_Radius) - stats(1).LineTraceTables{1,slice_to_plot}.Idealized_Radius(1)));
    ylim([stats(1).LineTraceTables{1,slice_to_plot}.Idealized_Radius(1)-yband-.5 stats(1).LineTraceTables{1,slice_to_plot}.Idealized_Radius(1)+yband+.5])
    pbaspect([4 1 1])
    xlim([ceil(min(stats(1).LineTraceTables{1,slice_to_plot}.Perimeter_Trace)) ceil(max(stats(1).LineTraceTables{1,slice_to_plot}.Perimeter_Trace))])
    xlabel('Perimeter (µm)')
    ylabel('Particle Radius (µm)')
    legend off
%     legend({'Fitted radius','Measured radius','Contact area'})
    
end