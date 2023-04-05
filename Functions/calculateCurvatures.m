function [edgecoor,edgecoor_sph,mc,gc] = calculateCurvatures(MPStats)
% Function to calculate mean and gaussian curvature over a rendered
% microparticle surface. Adapted from Vorselen et al 2020.
% Alex Settle & Miguel de Jesus
% Memorial Sloan Kettering Cancer Center
% Morgan Huse Laboratory, Department of Immunology
% 2023

[edgecoor,~,edgecoor_sph,TRI_nosph] = check_processing_status(MPStats(1));

[ThetaPhiRsmooth] = Smooth_spherical_mesh(edgecoor_sph,sqrt(1^2/pi),'arc',1,'mean');
[X,Y,Z] = sph2cart(ThetaPhiRsmooth(:,1),ThetaPhiRsmooth(:,2),ThetaPhiRsmooth(:,3));    

FV.faces    = MPStats(1).TRI_Connectivity_SPHbound;
FV.vertices = [X,Y,Z];
pc = GetCurvatures(FV,0);
gc = pc(1,:).*pc(2,:);
mc = mean(pc);


end