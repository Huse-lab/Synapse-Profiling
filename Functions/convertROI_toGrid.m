function [mc_i,gc_i,stain_i,scale] = convertROI_toGrid(MPStats,method,roi_data,grid_N)
% Convert MPStats data into a standardized grid in theta/phi
% space. 
%
% Will generate mc, gc and stain by default, unless it can't find stain
% data in which case it will just output mc,gc.
% 
% Alex Settle & Miguel de Jesus
% Memorial Sloan Kettering Cancer Center
% Morgan Huse Laboratory, Department of Immunology
% 2023
% Arguments
%  - MPStats(required): 
%  - method(required): determines whether contact area will be determined
%  by manual ROI or by previously-determined isincontact function in
%  MPStats. We recommend ROIs for most cases especilly when stain is dim.
%  - roi_data: if using ROI method, the roi_data structure
%  - grid_N: determine the resolution of gridded data, if not specified
%           will default to 50.
%
% Outputs
%   - mc_i = an NxN matrix representation of ROI containing mean curvature values
%   - gc_i = an NxN matrix representation of ROI containing gaussian curvature values
%   - stain_i = an NxN matrix representation of ROI containing intensity
%   values of the cellular stain (not particle)
%   - scale = the estimated scale of the grid in cartesian space, in µm

%Parse input arguments


%Get mc, gc and stain info from MPStats and ROI
iMP=1;
[edgecoor,centroidcoor,edgecoor_sph,TRI_nosph] = check_processing_status(MPStats(iMP));

if strcmp(method,'isincontact')
    contact_map = MPStats(iMP).isincontact;
elseif strcmp(method,'ROI')
    theta = edgecoor_sph(:,1); phi = edgecoor_sph(:,2);
    x = roi_data(iMP).position(:,1); y = roi_data(iMP).position(:,2); 
    contact_map = inpolygon(theta,phi,x,y); 
else
    fprintf('Error, method for selecting contact area not recognized')
    return
end

ThetaPhiR = edgecoor_sph;
R = ThetaPhiR(:,3);
R_contact = R(contact_map);

if isfield(MPStats,'stain_int')
    valid_stain = true;
    if strcmp(MPStats(iMP).Gridding,'equi')  
        stainint    = MPStats(iMP).stain_int;
        stainintint = MPStats(iMP).stain_int_integrated;
    elseif size(MPStats(iMP).stain_coor_cart,1)==size(MPStats(iMP).edgecoor_cart,1)
        stainint    = MPStats(iMP).stain_int;
        stainintint = MPStats(iMP).stain_int_integrated;
    elseif isfield(MPStats,'stain_int_irr')
        stainint    = MPStats(iMP).stain_int_irr;
        stainintint = MPStats(iMP).stain_int_integrated_irr;
    else
        try
            stainint    = MPStats(iMP).stain_int;
            stainintint = MPStats(iMP).stain_int_integrated;
            [~,stainint,~,stainintint] = ...
                Determine_Coverage(imgaussfilt3(MPStats(iMP).IMstain,1), MPStats(iMP).Weighted_Centroid'/MPStats(iMP).PixelSizeXY,...
                MPStats(iMP).EquivDiameter,MPStats(iMP).edgecoor_sph_wcatorigin(:,1),MPStats(iMP).edgecoor_sph_wcatorigin(:,2),[MPStats(iMP).PixelSizeXY MPStats(iMP).PixelSizeZ],...
                'R',MPStats(iMP).edgecoor_sph_wcatorigin(:,3)/MPStats(iMP).PixelSizeXY,'Gaussiansignal',MPStats(iMP).IMstain_isgaussian,'radialinterval' ,...
                MPStats(iMP).IMstain_radialinterval);
        catch
            warning(['Unable to get stain data for... ' MPStats(1).FileName])
            valid_stain = false;
        end    
    end
else
    warning(['MPStats does not contain stain_int ' MPStats(1).FileName])
    valid_stain = false;
end

[ThetaPhiRsmooth] = Smooth_spherical_mesh(edgecoor_sph,sqrt(1^2/pi),'arc',1,'mean');
[X,Y,Z] = sph2cart(ThetaPhiRsmooth(:,1),ThetaPhiRsmooth(:,2),ThetaPhiRsmooth(:,3));    

FV.faces    = MPStats(iMP).TRI_Connectivity_SPHbound;
FV.vertices = [X,Y,Z];
pc = GetCurvatures(FV,0);
gc = pc(1,:).*pc(2,:);
mc = mean(pc);

mc_contact = mc(contact_map);
gc_contact = gc(contact_map);
if(valid_stain)
    stain_int_contact = stainint(contact_map);
end

edgecoor_sph_contact = edgecoor_sph(contact_map,:);
edgecoor_contact = edgecoor(contact_map,:);

[~,~,mc_i] = interpolateCurvature(edgecoor_sph_contact(:,1), edgecoor_sph_contact(:,2),mc_contact,grid_N);   
mc_i(isnan(mc_i)) = 0;
    
[~,~,gc_i] = interpolateCurvature(edgecoor_sph_contact(:,1), ...
    edgecoor_sph_contact(:,2),gc_contact,grid_N);
gc_i(isnan(gc_i)) = 0;

if(valid_stain)
    [~,~,int_i] = ... 
            interpolateCurvature(edgecoor_sph_contact(:,1), ...
            edgecoor_sph_contact(:,2), double(stain_int_contact),grid_N);
    int_i(isnan(int_i)) = 0;
else
    int_i = [];
end

stain_i = int_i;
dist_matrix = pdist(edgecoor_contact);
scale = max(dist_matrix); % maximum diameter of contact area in µm



%% Functions 
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


end