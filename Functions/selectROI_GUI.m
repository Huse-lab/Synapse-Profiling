function ROI_data = selectROI_GUI(MPStats, display_method,draw_method)
% GUI to display worldmap view of particle and allow manual selection of
% the region of interest in theta/phi space
% 
% Alex Settle + Miguel de Jesus
% Memorial Sloan Kettering Cancer Center
% Morgan Huse Laboratory, Department of Immunology
% 2023
% 
% Arguments
% -MPStats:  An MPStats or MPRender structure output from ___
% -display_method: choose the values used for color scaling on the 
%            displayed particle
%         - 'Displacement'  - Edge Coordinate distance to particle centroid
%         - 'Stain' - Intensity of stain channel
% -draw_method: 'freehand' (hold click) vs. 'click' (click points)
%
% Output:
% -ROI_data: a matlab structure with the ROI label and polygon positions in
%   theta/phi space on the world map projection



% Parse arguments and determine options
display_options={'Displacement' 'Stain'};
color = find(strcmp(display_method,display_options));
if isempty(color)
    msg = "Error: invalid  display_method, choose 'Displacement'" + ...
        " or 'Stain'";
    error(msg)
    return
end

draw_options={'freehand' 'click'};
draw_type = find(strcmp(draw_method,draw_options));
if isempty(draw_type)
    msg = "Error: invalid  draw_method, choose 'freehand'" + ...
        " or 'click'";
    error(msg)
    return
end


%Load coordinate and stain information 
try 
    [edgecoor,centroidcoor,edgecoor_sph,TRI_nosph] = check_processing_status(MPStats(1));
catch
    error('Error: unable to load MPStats, check that data type is correct')
end


done = false;
while done == false
   % Plot the 2D surface
    fprintf('Drawing ROI for %s \n',MPStats(1).FileName);

    figure('Name',MPStats(1).FileName,'position',[0,300,2000,700])
    if color == 1
       R = edgecoor_sph(:,3);
       colormap(brewermap([],'RdBu')); 
       method = 'R';
    elseif color == 2
       R = MPStats(1).stain_int;
       stainmap = ones(256,3)*.05;
       stainmap(:,2) = linspace(.05,.9,256);
       colormap(stainmap)
       method = 'stain'; 
       try
           stainname = 'stain';
           stainint = MPStats(1).([stainname '_int']);
       catch
           error('Error: Stain values not present in MPStats structure')
       end
    else
       fprintf("Error. ")
       break 
    end
    ThetaPhiR = edgecoor_sph;
    ph = patch('faces',TRI_nosph,'vertices',[ThetaPhiR(:,1:2) R],...
        'facevertexcdata',R,'facecolor','flat','edgecolor','none');

   
if color==1
   caxis([prctile(R,.1) prctile(R,99.9)])
elseif color == 2
    caxis([prctile(stainint,1) prctile(stainint,99)]) 
end
   Set_2D_plot_opts
   add_colorscalebar()
              
   %Prompt user to select ROI
   if draw_type == 1
       roi = drawfreehand('FaceAlpha',0.2,'Color',[0.75 0.75 0.75]);
   elseif draw_type == 2
       roi = drawpolygon('FaceAlpha',0.2,'Color',[0.75 0.75 0.75]);
   end

   accept = questdlg('Are you satisfied with this ROI selection?', MPStats(1).FileName,'Yes','No, Try Again','Yes');
   % In case the user makes an error, they are given the option to redraw
   if strcmp(accept,'Yes')
       done = true;
   else
       close
   end

end
ROI_data = struct('position',{},'label',{});
ROI_data(1).position = roi.Position;
ROI_data(1).label = MPStats(1).FileName;
close


end 

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