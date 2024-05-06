function[Stats] = Update_Mask_GUI_MDJ(Stats,signals) 

% Some comment:
% Currently, you can use this function with only the deformation signal of
% the particle, BUT in this case snaking of the mask doesn't work.

if nargin == 1
    signals = '';
end

if isempty(signals) 
    if isfield(Stats,'stain_int')
        signals = {'stain_int'};
    end
    
    fnames = fieldnames(Stats);
    additional_stainchannel = contains(fnames,'stain_ch');
    if any(additional_stainchannel)
        fnames_addch = fnames(additional_stainchannel);
        ch_nmbr      = regexp(fnames_addch(2),'\d','match');
        signals{2} = ['stain_ch' ch_nmbr{1}{1} '_int'];
    end    
end

if isfield(Stats,'edgecoor_sph_aligned')
    edgecoor_sph     = Stats.edgecoor_sph_aligned;
    TRI_nosph        = Stats.TRI_Connectivity_aligned;
    if isempty(edgecoor_sph)
        edgecoor_sph = Stats.edgecoor_sph_wcatorigin;
        TRI_nosph    = Stats.TRI_Connectivity_wcatorigin;
    end
    shape_fieldname  = {'edgecoor_sph_aligned'};
    aligned_at_start = 1;
else
    edgecoor_sph     = Stats.edgecoor_sph_wcatorigin;
    TRI_nosph        = Stats.TRI_Connectivity_wcatorigin;
    shape_fieldname  = {'edgecoor_sph_wcatorigin'};
    aligned_at_start = 0;
end

if length(signals) < 2
    signals = [signals shape_fieldname];
end

% Set up plots
maskmap  = repmat(linspace(.05,.99,256)',1,3);
startmapFL = ones(256,3)*.05;
if contains(signals{1},'stain')
    stainmap  = startmapFL;
    stainmap(:,2)      = linspace(.05,.9,256);
else
    stainmap = colormap(brewermap([],'RdBu'));
end
if length(signals)>1
    if contains(signals{2},'stain')    
        stainmap2 = startmapFL;
        stainmap2(:,[1 3]) = [stainmap(:,2) stainmap(:,2)];
    else
        stainmap2 = colormap(brewermap([],'RdBu'));
    end
end
sh1 = []; sh2 = []; sh3 = [];

fh = figure;
fh.Position = [1580 312 762 888];

if isfield(Stats,'isincontact')
    mask             = Stats.isincontact;
    if ~aligned_at_start
        ReBase();
    end      
else
    mask = false(length(edgecoor_sph),1);
    if ~aligned_at_start
        warning('Starting without initial alignment');
    end
end
maskOG           = mask;
Snake_smoothing  = 1;

Make_Plots();

warning off backtrace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ROIh = [];
ROIcolor = 'r';

% Find components longitude, latitude and stress 
lons   = edgecoor_sph(:,1); 
lats   = edgecoor_sph(:,2);

% Add a button to allow switching between adding to and removing from mask the mask
bgh   = uibuttongroup('Units','normalized','Position',[0.8 0.65 0.15 0.1],'SelectionChangedFcn',@ Change_ROI_Color);
rembh = uicontrol(bgh,'Style','radiobutton','String','Remove','Units','normalized','Position',[0.1 0.55 0.8 0.35]);
addbh = uicontrol(bgh,'Style','radiobutton','String','Add','Units','normalized','Position',[0.1 0.1 0.8 0.35]);
if ~any(mask)
    addbh.Value = 1;
    ROIcolor = 'g';  
end

% Add a button to allow switching between adding to and removing from mask the mask
bgh2       = uibuttongroup('Units','normalized','Position',[0.8 0.55 0.15 0.1],'SelectionChangedFcn',@ Switch_ROI_Type);
circlebh   = uicontrol(bgh2,'Style','radiobutton','String','Circle' ,'Units','normalized','Position',[0.1 0.1 0.8 0.35]);
polygonbh  = uicontrol(bgh2,'Style','radiobutton','String','Polygon','Units','normalized','Position',[0.1 0.55 0.8 0.35]);

% Add a button to that performs realignment 
uicontrol('Parent',fh,'Style','pushbutton','String','Realign','Callback',@ReBase,'Units','normalized','Position',[0.8 0.5 0.15 0.05]);

if length(signals) > 1
    % And a button that performs the snake algorithm
    uicontrol('Parent',fh,'Style','text','String','Enter desired smoothing strength','Units','normalized','Position',[0.8 0.38 0.15 0.05]);
    tiph = uicontrol('Parent',fh,'Style','edit','String','3','Callback',@Update_smoothing,'Units','normalized','Position',[0.8 0.35 0.15 0.05]);
    uicontrol('Parent',fh,'Style','pushbutton','String','Snake my Mask!','Callback',@Snake_this_Mask,'Units','normalized','Position',[0.8 0.3 0.15 0.05]);
end

% Add a button that removes the entire mask or resets it to the original
uicontrol('Parent',fh,'Style','pushbutton','String','Clear Mask(!)','Callback',@Clear_Mask ,'Units','normalized','Position',[0.8 0.175 0.15 0.05],'BackgroundColor',[0.95 0.3 0.3]);
uicontrol('Parent',fh,'Style','pushbutton','String','Invert Mask'  ,'Callback',@Invert_Mask,'Units','normalized','Position',[0.8 0.125 0.15 0.05]);
uicontrol('Parent',fh,'Style','pushbutton','String','Reset Mask'   ,'Callback',@Reset_Mask ,'Units','normalized','Position',[0.8 0.075 0.15 0.05]);
uicontrol('Parent',fh,'Style','pushbutton','String','Done!'        ,'Callback',@Close_Fig  ,'Units','normalized','Position',[0.8 0.025 0.15 0.05],'BackgroundColor',[0.3 0.95 0.3]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wait until the user closes the figure to finalize

waitfor(fh);

if ~isfield(Stats,'isincontact')
    latestmaskchangesapplied = 0;
elseif ~all(Stats.isincontact == mask)
    latestmaskchangesapplied = 0;
else
    latestmaskchangesapplied = 1;
end

if ~latestmaskchangesapplied
    applylastchanges = questdlg(['Latest changes made in plot not yet applied to mask. '...
        'Do you want to apply changes now?'],'Apply latest changes?','Yes','No','Yes');
    if applylastchanges
        ReBase;
        close(gcf);
    end
end


if nargout == 0 && ~all(maskOG == mask)
    warning('No direct function output. Stats available as "ans".')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callback Functions

    % Function that allows adding or removing mask segments
    function [] = Add_or_Remove_from_Mask(~,~)
               
        % Make a ROI
        if circlebh.Value
            ROIh = drawcircle(gca,'Color',ROIcolor,'LineWidth',1,'DrawingArea','unlimited');
        elseif polygonbh.Value
            ROIh = drawpolygon(gca,'Color',ROIcolor,'LineWidth',1,'DrawingArea','unlimited');
        end
        
        % Prevent user from initiating multiple ROIs
         Allow_new_ROI(0);
        
        % Wait until the user is finished
        try
            disp('Once you have finalized your mask section, press enter to apply')
            w = waitforbuttonpress;
            % Wait for keyboard press, ignoring mouseclicks
            while w == 0;     w = waitforbuttonpress;        end 
        catch
            return;
        end
        
        % Find coordinates in circle
        if circlebh.Value && strcmpi(ROIh.Type,'images.roi.Circle')
            inROI = sqrt((lons-ROIh.Center(1)).^2+(lats-ROIh.Center(2)).^2)<ROIh.Radius;  
        elseif polygonbh.Value && strcmpi(ROIh.Type,'images.roi.polygon')
            vertexlist = ROIh.Position;
            inROI = inpolygon(lons,lats,vertexlist(:,1),vertexlist(:,2));
        else
            delete(ROIh)
            return;
        end
             
        % Add or remove area from mask
        if addbh.Value           
            mask(inROI) = 1;            
        elseif rembh.Value
            mask(inROI) = 0;
        end       
        
        % Remove circle and update plot
        delete(ROIh)
        Make_Plots(1);
        
    end


    % Switch between adding circular or polygonal sections to the mask 
    function [] = Switch_ROI_Type(~,~)
        
        delete(ROIh)
        Allow_new_ROI(1);
        
    end


    % Change between green for adding and red for removing sections from the mask 
    function [] = Change_ROI_Color(~,~)
        
        if addbh.Value           
            ROIcolor = 'g';   
        elseif rembh.Value
            ROIcolor = 'r';
        end     
        
        if ~isempty(ROIh)
            if isvalid(ROIh)
                ROIh.Color = ROIcolor; 
            end
        end
        
    end


    % Here we realign the particle and incorporate the changes in the
    % output structure
    function [] = ReBase(~,~)
       
       % Determine the position of the base and align cups 
       Stats.isincontact = mask;    
       if ~all(mask(:) == mask(1))
            Stats = Determine_base_position_and_align_MDJ(Stats,'BaseLat',-pi/2,'BaseColongitude',0);
            edgecoor_sph     = Stats.edgecoor_sph_aligned;                  
            lons             = edgecoor_sph(:,1); 
            lats             = edgecoor_sph(:,2);
            Stats.Fraction_engulfed = Calc_Frac_engulfed(lons,lats,mask,length(mask));
            TRI_nosph        = Stats.TRI_Connectivity_aligned;
            Make_Plots;
       else
           warning('Current mask does not distinguish any foreground and background. Rebasing not possible')
       end
        
    end        


    % Try optimizing the user defined mask using localized segmentation
    function[] = Snake_this_Mask(~,~)
        
        if all(mask(:) == mask(1))
            warning(['Current mask implies everything is foreground or everything is ',...
                'background. Make sure to add or remove some fragments to the mask'])
            return;
        end       
        
        % Store the old mask in case snaking gives unsatisfactory results
        oldmask = mask;
              
        % Obtain edge coordinates
        theta_irr = Stats.edgecoor_sph_wcatorigin(:,1);
        phi_irr   = Stats.edgecoor_sph_wcatorigin(:,2);
        
        % Use a smoothed version of the image for snaking
        regstainpresent = 0;
        if isfield(Stats,'edgecoor_gridsize')
            if ~isempty(Stats.edgecoor_gridsize)
                gridsize1 = Stats.edgecoor_gridsize(1);
                gridsize2 = Stats.edgecoor_gridsize(2);
                if isfield(Stats,'IMstain_reg')
                    pixintgrid = reshape(Stats.IMstain_reg,gridsize1,gridsize2);
                elseif Stats.stain_indicates_contact
                    pixintgrid = reshape(Stats.IMstain2_reg,gridsize1,gridsize2);
                else
                    pixintgrid = reshape(Stats.IMstain1_reg,gridsize1,gridsize2);
                end
                
                Theta = repmat(linspace(-pi,pi * ( 1 - (2 / gridsize1 ) ),gridsize1),1,gridsize2);
                Phi   = reshape(repmat(linspace(-.5*pi,.5*pi,gridsize2),gridsize1,1),1,[]); 
                
                regstainpresent = 1;
            end
        end
        
        if ~regstainpresent
            [Theta,Phi,gridsize] = Pack_regular_Points_on_Sphere(round(sqrt(length(Stats.edgecoor_cart_wcatorigin)*2)));
            gridsize1 = gridsize(1);
            gridsize2 = gridsize(2);
            pixints = Interpolate_spherical_surface_MDJ(Theta,Phi,theta_irr,phi_irr,Stats.stain_int);
            pixintgrid = reshape(pixints,gridsize1,gridsize2);
        end
        
        if Snake_smoothing > 0;     pixintgrid = imgaussfilt(pixintgrid,Snake_smoothing);       end  

        % interpolate the most recent mask on a sphere    
        IPmask = Interpolate_spherical_surface_MDJ(Theta,Phi,theta_irr,phi_irr,mask) > 0.5;
        IPmask = reshape(IPmask,gridsize1,gridsize2);
        
        % Apply localized segmentation (the snake)
        [newmask,phi] = localized_seg(log(pixintgrid),IPmask,200,5,.1,[],[],1);
    
        % Obtain a mask corresponding to our equidistant points
        GI = griddedInterpolant({[Theta(1:gridsize1) pi],Phi(1:gridsize1:end)},double([newmask; newmask(1,:)]));
        iscovered = logical(round(GI(theta_irr,phi_irr)));  
        iscovered(iscovered<0) = 0;
        mask = iscovered;
        Make_Plots(1);
        
        % Ask the user if the snake yielded satisfactory results
        answer = questdlg('Do you want to apply or undo snaking?','Satisfactory snaking?','Apply','Undo','Undo');
        
        % Apply or discard changes
        if strcmpi(answer,'Apply')               
                Stats.actcont_phi = phi;
                Stats.actcont_mask  = single(mask);
                ReBase;
                Make_Plots();            
        elseif strcmpi(answer,'Undo')                
                mask = oldmask;
                Make_Plots(1);                
        end                   
    end


    % Update the smoothing parameter for the localized segmentation algorithm
    function [] = Update_smoothing(~,~)
        
        Snake_smoothing = str2double(tiph.String);
        
    end


    % Clear mask completely, allows drawing a mask from scratch
    function[] = Clear_Mask(~,~)
        
        mask = false(size(mask));
        Make_Plots(1);
        addbh.Value = 1;
        ROIcolor = 'g';      
        
    end


    % Invert mask
    function[] = Invert_Mask(~,~)
        
        mask = ~mask;
        Make_Plots(1);
        
    end


    % Undo all user based changes
    function[] = Reset_Mask(~,~)
        
        mask = maskOG;
        Make_Plots(1);
        
    end


    function[] = Close_Fig(~,~)
        
        delete(fh);
        
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other functions

    % Make Plots
    function[] = Make_Plots(varargin)
        
        if nargin == 0
            whichplots = [1 (1:length(signals))+1];
        else
            whichplots = [varargin{:}];
        end
        
        if any(whichplots==1)
            sh1 = subplot(3,1,1);
            delete(sh1.Children);
            ph1 = patch('faces',TRI_nosph,'vertices',edgecoor_sph,'facevertexcdata',mean(mask(TRI_nosph),2),'facecolor','flat','edgecolor','none');
            colormap(maskmap)
            caxis([0 1])
            Set_2D_plot_opts()
            title('Mask')
            ph1.HitTest = 'off';
            sh1.PickableParts = 'all';
        end

        if any(whichplots==2)
            sh2 = subplot(3,1,2);
            delete(sh2.Children);
            stainint = Stats.(signals{1});
            % This is a quick fix for backward compatibility with old
            % MPStats. It should go in a future version.
            if length(stainint)~=length(edgecoor_sph)
                [~,stainint] = Determine_Coverage(imgaussfilt3(Stats.IMstain,1), Stats.Weighted_Centroid'/Stats.PixelSizeXY,...
                    Stats.EquivDiameter,Stats.edgecoor_sph_wcatorigin(:,1),Stats.edgecoor_sph_wcatorigin(:,2),[Stats.PixelSizeXY Stats.PixelSizeZ],...
                    'R',Stats.edgecoor_sph_wcatorigin(:,3)/Stats.PixelSizeXY,'Gaussiansignal',Stats.IMstain_isgaussian,'radialinterval',Stats.IMstain_radialinterval);
                stainint = stainint{1};
            end
            if size(stainint,2) == 3;       stainint = stainint(:,3);       end
            ph2 = patch('faces',TRI_nosph,'vertices',edgecoor_sph,'facevertexcdata',mean(stainint(TRI_nosph),2),'facecolor','flat','edgecolor','none');
            colormap(stainmap)
            caxis([prctile(stainint,1) prctile(stainint,99)])
            Set_2D_plot_opts   
            title(signals{1},'interpreter','none')
            ph2.HitTest = 'off';
            sh2.PickableParts = 'all';
        end

        if any(whichplots==3)
            sh3 = subplot(3,1,3);
            delete(sh3.Children);
            stainint2 = Stats.(signals{2});
            if size(stainint2,2) == 3;      stainint2 = stainint2(:,3);     end
            ph3 = patch('faces',TRI_nosph,'vertices',edgecoor_sph,'facevertexcdata',mean(stainint2(TRI_nosph),2),'facecolor','flat','edgecolor','none');
            colormap(stainmap2)
            caxis([prctile(stainint2,1) prctile(stainint2,99)])
            Set_2D_plot_opts    
            title(signals{2},'interpreter','none')
            ph3.HitTest = 'off';
            sh3.PickableParts = 'all';
        end    
        
        % Set windowkeypressfcn for fh
        Allow_new_ROI(1)     
        
    end


    % Set keypressfunctions for figure
    function[] = Allow_new_ROI(allowed)
        
        if allowed
            sh1.ButtonDownFcn = @(x,event) Add_or_Remove_from_Mask(x,event);
            sh2.ButtonDownFcn = @(x,event) Add_or_Remove_from_Mask(x,event);
            sh3.ButtonDownFcn = @(x,event) Add_or_Remove_from_Mask(x,event);
        else
            sh1.ButtonDownFcn = [];
            sh2.ButtonDownFcn = [];
            sh3.ButtonDownFcn = [];
        end        
    end


    % Set 2D plot options
    function[] = Set_2D_plot_opts()

        view(3);
        axis off
        axis equal
        view(0,90);
        set(gca,'Clipping','off')
        %freezeColors;
        axtoolbar();

    end

end