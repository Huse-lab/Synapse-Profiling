function [Stats] = Determine_base_position_and_align(Stats,varargin)

% Determine the position of the base of the cup and align cups
ThetaPhi  = cellfun(@(x) x(:,[1 2]),{Stats.edgecoor_sph_wcatorigin},'UniformOutput',false);

% The centroid of the contact area is taken
[~,~,~,theta_base,phi_base] = Determine_Centroid(ThetaPhi,{Stats.isincontact});

[ThetaRot,PhiRot] = Align_spheres_by_Rotation(ThetaPhi,double([theta_base phi_base]));
% Note that it is still possible to rotate around the base - opposite apex axis

NMPs = length(Stats);

% Redo triangulation
for iMP = 1:NMPs 

    % Determine theta, phi, R
    theta     = wrapToPi(ThetaRot{iMP}-.5*pi); % We rotate theta away from the center for visualization on 2D maps
    phi       = PhiRot{iMP};
    r         = Stats(iMP).edgecoor_sph_wcatorigin(:,3);
       
    Stats(iMP).edgecoor_sph_aligned =     [theta,phi,r];
    Stats(iMP).TRI_Connectivity_aligned = delaunay(theta,phi);
    
    [x,y,z] = sph2cart(theta,phi,r);
    Stats(iMP).edgecoor_cart_aligned = [x,y,z];
    
    % Also save the location of the cupbase
    [~,minloc] = min(pdist2([theta,phi],[.5*pi,0]));
    r_base = r(minloc);
    Stats(iMP).CupBase = [theta_base(iMP,:) phi_base(iMP,:) r_base];
    
end

Stats = orderfields(Stats);