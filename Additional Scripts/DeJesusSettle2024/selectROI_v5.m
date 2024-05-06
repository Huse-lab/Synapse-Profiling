% ROI selection script on MPRender to select 

% This code loops through a selected data folder, opens the MPStats file
% for each selected frame, and allows the user to defined the region of
% interest (ROI) by drawing an arbitrary polygon. It will create a .mat
% file with the select points on the polygon which can be used later for
% analysis.

%updated 1/8/2021 to v3 to allow option to show radius or stain when
%selecting ROIs
% Updated by MDJ 2021-01-08 to allow free hand drawing, for ease of use with a pen
% tablet.
% Updated by MDJ 2021-01-11 to adjust file name without strrep, because
% "MPRender" appears in the folder name as well.
% Updated by MDJ 2021-03-09 to uniformly display the maps for selecting
% "view(0,90)"; also added quality of life user input lines to ask about R
% map vs. stain map.
% Updated 2021-07-15, MDJ: To adjust to new Daan code.

%% Select Folder of interest

file_path = uigetdir;
file_list = dir(file_path);

%Determine what color will be displayed, 0 = displacement, 1 = stain
color = input("Visualize particles using displacement or cell stain? (0: R map, 1: stain map) "); 

%% Loop through MPStats files

warning('off','all')

for i = 1:length(file_list)
   if isempty(strfind(file_list(i).name,'MPRender'))
       continue
   end
   disp(file_list(i).name)  
   file_name = strcat(file_path, '/', file_list(i).name);
   % changed how filenames are defined, to avoid strrepping the folder
   % names: replace only the last instance of "MPRender_pool" with
   % "ROI_pool"
   index_match = strfind(file_name,'MPRender');
   match_2 = index_match(end);
   roi_filename = [file_name(1:match_2-1),'ROI_pool',file_name(match_2+length('MPRender_pool'):end)];
   %roi_filename = strrep(file_name,'MPRender','ROI');
   
   %check if ROI file already exists
   if isfile(roi_filename)
       overwrite = questdlg('An ROI.mat file already exists for this image. Would you like to overwrite?', ... 
           file_list(i).name,'Yes','No','No');
       if strcmp(overwrite,'No')
           continue
       end
   end
   load(file_name);
   roi_position = struct('label',{},'position',{});
   
   
   if ~exist('iMP','var'); iMP = 1; end
   for k = 1:length(MPStats)
       [edgecoor,centroidcoor,edgecoor_sph,TRI_nosph] = check_processing_status(MPStats(k));
       
       done = false;
       while done == false
           % Plot the 2D surface
           fprintf('Drawing ROI for %s \n',MPStats(k).FileName);

           figure('Name',MPStats(k).FileName,'position',[0,300,2000,700])
           if color == 0
               R = edgecoor_sph(:,3);
               colormap(brewermap([],'RdBu')); 
               method = 'R';
           elseif color == 1
               R = MPStats(k).stain_int;
               stainmap = ones(256,3)*.05;
               stainmap(:,2) = linspace(.05,.9,256);
               colormap(stainmap)
               method = 'stain';                           
           else
               fprintf("Error. ")
               break 
           end
           
                ThetaPhiR = edgecoor_sph;
                stainname = 'stain';
                stainint = MPStats(k).([stainname '_int']);
            ph = patch('faces',TRI_nosph,'vertices',[ThetaPhiR(:,1:2) R],'facevertexcdata',R,'facecolor','flat','edgecolor','none');
%               'R' instead of 'stainint' if doing selections on R, vice
%               versa if actin stain.
            view(3)
            axis off
            axis equal
            view(0,90);
            set(gca,'Clipping','off')
           
        if color==0
%            h = trisurf(triangulation(TRI_nosph,[edgecoor_sph(:,[1 2]),R]),R,'LineStyle','none');
           caxis([prctile(R,.1) prctile(R,99.9)])
%             colormap(stainmap)
        else
            caxis([prctile(stainint,1) prctile(stainint,99)]) 
        end
           Set_2D_plot_opts
           add_colorscalebar()
           
           % added this 2021-03-09 for dataset IMG117
           view(0,90)
                      
           %Prompt user to select ROI
           roi = drawfreehand('FaceAlpha',0.1,'Color',[0.75 0.75 0.75]);        
           accept = questdlg('Are you satisfied with this ROI selection?', file_list(i).name,'Yes','No, Try Again','Yes');
           if strcmp(accept,'Yes')
               done = true;
           else
               close
           end
       
       end
       
       roi_position(k).position = roi.Position;
       roi_position(k).label = MPStats(k).FileName;
       close
   end
   
   if ~contains(roi_filename,'ROI') % doublecheck, only overwriting ROIs
       printf('Error, trying to overwrite wrong file')
   else
       fprintf(strcat('Saving','->',roi_filename,'\n'))
       save(roi_filename,'roi_position')
   end
   
   close
   
       
end
 
disp('Done drawing ROIs for this folder.')

%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select particle for plot
function iMP = Select_particle_for_plot(MPStats,iMP)

    NMPs = length(MPStats);

    if NMPs > 1

        iMP = str2double(inputdlg('which particle do you want to visualize?',...
            'Select particle',[1 50],{num2str(iMP)}));

        if iMP > NMPs
            error(['Particle ' num2str(iMP) ' doesn''t exist']);
        end

    else

        iMP = 1;

    end
end

% Set 3D plots options
function[] = Set_3D_plot_opts()

    zoom(1.25)
    axis equal
    axis off
    set(gca,'CameraViewAngleMode','Manual','Clipping','off')
    view(0,0)

end

% Set 2D plot options
function[] = Set_2D_plot_opts()

    axis off
    axis equal
    xlabel('Theta')
    ylabel('Phi')
    grid off
    view(0,-90);
    set(gca,'Clipping','off')

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

function drawcircle(x,y,r)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp,'g-','LineWidth',3);
end


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