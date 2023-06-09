function LPA_Stats3D = calculateBulkDeformation(MPStats,roi_data,num_slices)
% Slice-based analysis of 3D-deformation of particle induced by cell
% 
%
% Alex Settle & Miguel de Jesus
% Memorial Sloan Kettering Cancer Center
% Morgan Huse Laboratory, Department of Immunology
% 2023
%
% Arguments:
%      -MPStats: An MPStats or MPRender structure 
%      -roi_data: corresponding ROI data, generated by selectROI_GUI
% Output: 
%      -LPA_Stats3D: A Data structure containing deformation calculation results     
%           edgecoor_cart_aligned
%                         FileName: str, the corresponding image being analyzed
%            edgecoor_cart_aligned: Nx3 double, cartesian coordinates of particle edge
%             edgecoor_sph_aligned: NX3 double, spherical coordiinates of particle edge
%                        stain_int: Nx1 array, Cellular stain intensity at particle edge
%                     roi_position: Nx2 double, bounding polygon defined by ROI 
%                     ROI_Z_bounds: double, maximum and minimum Z range of ROI
%                 ROI_Theta_Bounds: double, max and minimum theta coordinate for
%                                   ROI, used for plotting
%                      Num_ZSlices: integer, Number of subdivisions used to
%                                     calculate deformaton
%                  LineTraceTables: cell-array, perimeter around the
%                                    particle at each z-section calculated.
%                                    Useful for plotting.
%                                   Within each table:
%                                       -Theta - theta coordinate
%                                       -Perimeter_Trace - perimeter
%                                       coordinate around the particle
%                                       -Idealized_Radius - the estimated
%                                       radius of the particle at that
%                                       specific z-slice, estimated from
%                                       non-contacted portion of particle
%                                       -Measured_Radius - the measured
%                                       distance of particle edge to center
%                                       at each point on the
%                                       Perimeter/Theta continuum.
%                                   %to generate line profile:
%                                   plot Theta/Perimeter vs Measured_Radius  
%                   CenterDepth_um: scalar, the maximum indentation from an
%                                   idealized circle at the center of contact,
%                                   estimated from ± 15% from center
%       Center_Slice_Deviation_ROI: scalar, The area deviation at center of
%                                   contact, only within ROI
%     Center_Slice_Deviation_Total: the total area deviation at the center
%                                   of contact
%             Volume_Deviation_ROI: scalar, Volumetric Deformation within 
%                                   the ROI in µm^3
%           Volume_Deviation_Total: scalar, Total Volumetric Deformation,
%                                     µm^3
%               Est_Force_Hertz_pN: Estimated force of deformation based on
%                                    Hertz Model of contact (see methods),
%                                    in pico-newton
%                Est_BulkStress_Pa: Estimated Bulk stress based on dV, in
%                                   pico-newton
%                 Est_BulkForce_pN: Estimated total Bulk force based on dV,
%                                   in pico-newton
%                   ParticleRadius: Estimated unindented particle radius, µm
%            InitialParticleVolume: Estimated unindented particle Volume,
%                                   µm^3
%                         dVStrain: Unitless Strain, calculated as
%                                  volume deformation / total volume


i=1;

%Initialize LPA_Stats3D and get neceessary info from MPStats and ROI
LPA_Stats3D = struct('FileName',{},'edgecoor_cart_aligned',[],'edgecoor_sph_aligned',[],'stain_int',[],'roi_position',[],'ROI_Z_bounds',[],'ROI_Theta_Bounds',[],'Num_ZSlices',[],...
        'LineTraceTables',[],'CenterDepth_um',[],'Center_Slice_Deviation_ROI',[],'Center_Slice_Deviation_Total',[], ...
        'Volume_Deviation_ROI',[],'Volume_Deviation_Total',[],'Est_Force_Hertz_pN',[],...
        'Est_BulkStress_Pa',[],'Est_BulkForce_pN',[],'ParticleRadius',[],'InitialParticleVolume',[],...
        'dVStrain',[]);
LPA_Stats3D(i).FileName = MPStats(i).FileName;
LPA_Stats3D(i).edgecoor_cart_aligned = MPStats(i).edgecoor_cart_aligned;
LPA_Stats3D(i).edgecoor_sph_aligned = MPStats(i).edgecoor_sph_aligned;
LPA_Stats3D(i).stain_int = MPStats(i).stain_int;
LPA_Stats3D(i).roi_position = roi_data(i).position;
volumeStats = calculate_volumeStats(LPA_Stats3D(i),num_slices);
LPA_Stats3D(i) = volumeStats;





end


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
    
    volumeStats.CenterDepth_um = centerslice_Stats.MaximumDepth_um;
    volumeStats.Center_Slice_Deviation_ROI = centerslice_Stats.Integrated_Deviation_ROI;
    volumeStats.Center_Slice_Deviation_Total = centerslice_Stats.Integrated_Deviation_Total;
    

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
        

end