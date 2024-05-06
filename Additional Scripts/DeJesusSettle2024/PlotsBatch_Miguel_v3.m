%% Initialization

% Also plot curvatures and stain
plotcurv      = 1;
plotstain     = 1;
plotmask      = 0;

ntot = length(MPStats);
rfig = ceil(sqrt(ntot));
cfig = ceil(ntot/rfig);
subplotn = [rfig,cfig];

%% Plot triangulated surface of particles

if plotstain
    stainname = 'stain';
    fnames = fieldnames(MPStats);
    second_channel_present = any(contains(fnames,'stain_ch'));
    if second_channel_present
        plotstain2 = 1;
        fnames_secondchannel = fnames(contains(fnames,'stain_ch'));
        stainname2 = fnames_secondchannel{1}(strfind(fnames_secondchannel{1},'stain_ch')+(0:8));
    else
        plotstain2 = 0;
    end
end

nbatches = ceil(length(MPStats)/ntot);
AH = cell(nbatches,1);
if exist('plotstain','var') 
    if plotstain == 1
        AH2 = cell(nbatches,1);
        AH3 = cell(nbatches,1);
    end
end
if exist('plotcurv','var')    
    if plotcurv == 1
        AH4 = cell(nbatches,1);
        AH5 = cell(nbatches,1);
    end
end
if exist('plotmask','var')
    if plotmask == 1
        AH6 = cell(nbatches,1);
    end
end

for jbatch = 1:nbatches
    
    fh = figure('Name','Distance to Centroid Sphere');
    fh.Position = [1100 80 1064 1064];
    ah = tight_subplot(subplotn(1),subplotn(2),[.01 .01],[0.01 0.01],[0.01 0.01]);
    
    if exist('plotstain','var')
        if plotstain == 1
            if isfield(MPStats,[stainname '_int'])
                % Actin
                fh2 = figure('Name','LifeAct');
                fh2.Position = [1100 80 1064 1064];
                ah2 = tight_subplot(subplotn(1),subplotn(2),[.01 .01],[0.01 0.01],[0.01 0.01]);
                stainmap = ones(256,3)*.05;
                stainmap(:,2) = linspace(.05,.9,256);
            end

            if plotstain2
                if isfield(MPStats,[stainname2 '_int'])
                    % IO_Stain
                    fh3 = figure('Name','Stain');
                    fh3.Position = [1100 80 1064 1064];
                    ah3 = tight_subplot(subplotn(1),subplotn(2),[.01 .01],[0.01 0.01],[0.01 0.01]);

                    stainmap2 = ones(256,3)*.05;
                    stainmap2(:,1) = linspace(.05,.9,256);
                    stainmap2(:,3) = linspace(.05,.9,256);
                end
            end
        end
    end
    
    if exist('plotcurv','var')
        if plotcurv == 1
            
            fh4 = figure('Name','Mean Curvature deviation');
            fh4.Position = [1200 150 1064 1064];
            ah4 = tight_subplot(subplotn(1),subplotn(2),[.01 .01],[0.01 0.01],[0.01 0.01]);
            
            fh5 = figure('Name','Gaussian Curvature');
            fh5.Position = [1200 150 1064 1064];
            ah5 = tight_subplot(subplotn(1),subplotn(2),[.01 .01],[0.01 0.01],[0.01 0.01]);
            
        end
    end
    
    if exist('plotmask','var')
        if plotmask == 1           
            fh6 = figure('Name','Mask');
            fh6.Position = [1200 150 1064 1064];
            ah6 = tight_subplot(subplotn(1),subplotn(2),[.01 .01],[0.01 0.01],[0.01 0.01]);                       
        end
    end
    
    for iMP = (1+(jbatch-1)*ntot):min(ntot*jbatch,length(MPStats))
        [~,centroidcoor,edgecoor_sph,TRI_nosph] = check_processing_status(MPStats(iMP));
        ThetaPhiR = edgecoor_sph;
        R = ThetaPhiR(:,3);
        [X,Y,Z] = sph2cart(ThetaPhiR(:,1),ThetaPhiR(:,2),ThetaPhiR(:,3));

        % Plot the shape
        set(groot,'CurrentFigure',fh)
        fh.CurrentAxes = ah(iMP-(jbatch-1)*ntot); 
        Make_3D_plot(MPStats(iMP).TRI_Connectivity_SPHbound,[Z,Y,X],R);
        colormap(brewermap([],'RdBu'))
        caxis(median(R)+[-.5 .5])

        if exist('plotstain','var')
            if plotstain == 1
                %Plot the LifeAct
                if isfield(MPStats,[stainname '_int'])
                    set(groot,'CurrentFigure',fh2)
                    fh2.CurrentAxes = ah2(iMP-(jbatch-1)*ntot); 
                    stainint = MPStats(iMP).([stainname '_int']);
                    Make_3D_plot(MPStats(iMP).TRI_Connectivity_SPHbound,[Z,Y,X],stainint);
                    colormap(stainmap)
                    caxis([prctile(MPStats(iMP).stain_int,1) prctile(MPStats(iMP).stain_int,99.5)])
                end

                %Plot the stain
                if plotstain2
                    if isfield(MPStats,[stainname2 '_int'])
                        set(groot,'CurrentFigure',fh3)
                        fh3.CurrentAxes = ah3(iMP-(jbatch-1)*ntot); 
                        stainint = MPStats(iMP).([stainname2 '_int']);
                        Make_3D_plot(MPStats(iMP).TRI_Connectivity_SPHbound,[Z,Y,X],stainint);
                        colormap(stainmap2)
                        caxis([prctile(stainint,1) prctile(stainint,99)])
                    end
                end
            end
        end
        
        if exist('plotmask','var')
            if plotmask == 1
                %Plot the Mask
                set(groot,'CurrentFigure',fh6)
                fh6.CurrentAxes = ah6(iMP-(jbatch-1)*ntot); 
                stainint = MPStats(iMP).isincontact;
                Make_3D_plot(MPStats(iMP).TRI_Connectivity_SPHbound,[Z,Y,X],stainint);
                colormap(gray)
                caxis([0 1])
            end
        end
        
        if exist('plotcurv','var')    
            if plotcurv == 1
                
                [ThetaPhiRsmooth] = Smooth_spherical_mesh(edgecoor_sph,sqrt(1^2/pi),'arc',1,'mean');
                [X,Y,Z] = sph2cart(ThetaPhiRsmooth(:,1),ThetaPhiRsmooth(:,2),ThetaPhiRsmooth(:,3));    
            
                FV.faces    = MPStats(iMP).TRI_Connectivity_SPHbound;
                FV.vertices = [X,Y,Z];
                pc = GetCurvatures(FV,0);
                gc = pc(1,:).*pc(2,:);
                mc = mean(pc);

                set(groot,'CurrentFigure',fh4)
                fh4.CurrentAxes = ah4(iMP-(jbatch-1)*ntot); 
                Make_3D_plot(MPStats(iMP).TRI_Connectivity_SPHbound,[Z,Y,X], mc - median(mc));
                colormap(greenmagenta)
                caxis([-0.3 0.3])

                set(groot,'CurrentFigure',fh5)
                fh5.CurrentAxes = ah5(iMP-(jbatch-1)*ntot); 
                Make_3D_plot(MPStats(iMP).TRI_Connectivity_SPHbound,[Z,Y,X],gc);
                colormap(turquoisebrown)
                caxis([-0.18 0.18])
                
            end
        end
        
    end
    
    AH{jbatch}  = ah ;
    if exist('plotstain','var') 
        if plotstain == 1
            if isfield(MPStats,[stainname '_int'])
                AH2{jbatch} = ah2;
                if plotstain2
                    if isfield(MPStats,[stainname2 '_int'])
                        AH3{jbatch} = ah3;
                    end
                end
            end
        end
    end
    if exist('plotcurv','var')    
        if plotcurv == 1
            AH4{jbatch} = ah4;
            AH5{jbatch} = ah5;
        end
    end
    if exist('plotmask','var')
        if plotmask == 1
            AH6{jbatch} = ah6;
        end
    end
end

if rem(iMP,ntot)~=0
    for iax = ntot:-1:rem(iMP,ntot)+1
        delete(ah(iax))
        ah(iax) = [];
        if exist('plotstain','var')
            if plotstain == 1
                if isfield(MPStats,[stainname '_int']) ;    delete(ah2(iax));   ah2(iax) = [];       end
                if plotstain2
                    if isfield(MPStats,[stainname2 '_int']);    delete(ah3(iax));   ah3(iax) = [];       end
                end
            end
        end
        if exist('plotmask','var')
            if plotmask == 1;   delete(ah6(iax));   ah6(iax) = [];      end
        end
        if exist('plotcurv','var')    
            if plotcurv == 1   
                delete(ah4(iax));   delete(ah5(iax));
                ah4(iax) = []   ;   ah5(iax) = []   ;
            end
        end
    end    
end

AH{jbatch}  = ah ;
if exist('plotstain','var') 
    if plotstain == 1
        if isfield(MPStats,[stainname '_int'])
            AH2{jbatch} = ah2;
            if plotstain2
                if isfield(MPStats,[stainname2 '_int'])
                    AH3{jbatch} = ah3;
                end
            end
        end
    end
end
if exist('plotcurv','var')    
    if plotcurv == 1
        AH4{jbatch} = ah4;
        AH5{jbatch} = ah5;
    end
end
if exist('plotmask','var')
    if plotmask == 1
        AH6{jbatch} = ah6;
    end
end

%% Plot triangulated 2D surface of particles

% Use mollweid projections (if not, equirectangular map projections are used)
mollweid_proj = 0;

if exist('plotstain','var')   
    if plotstain
        stainname = 'stain';
        fnames = fieldnames(MPStats);
        second_channel_present = any(contains(fnames,'stain_ch'));
        if second_channel_present
            plotstain2 = 1;
            fnames_secondchannel = fnames(contains(fnames,'stain_ch'));
            stainname2 = fnames_secondchannel{1}(strfind(fnames_secondchannel{1},'stain_ch')+(0:8));
        else
            plotstain2 = 0;
        end
    end
end

nbatches = ceil(length(MPStats)/ntot);
AH = cell(nbatches,1);
if exist('plotstain','var') 
    if plotstain == 1
        AH2 = cell(nbatches,1);
        AH3 = cell(nbatches,1);
    end
end
if exist('plotcurv','var')    
    if plotcurv == 1
        AH4 = cell(nbatches,1);
        AH5 = cell(nbatches,1);
    end
end
if exist('plotmask','var')
    if plotmask == 1
        AH6 = cell(nbatches,1);
    end
end

for jbatch = 1:nbatches
    
    fh = figure('Name','Distance to Centroid Map');
    fh.Position = [1100 80 1064 1064];
    ah = tight_subplot(subplotn(1),subplotn(2),[.01 .01],[0.01 0.01],[0.01 0.01]);
        
    if exist('plotstain','var')    
        if plotstain
            if isfield(MPStats,[stainname '_int'])
                % Actin
                stainlabel = input("What is the stain on this cell? ",'s');
                fh2 = figure('Name',stainlabel);
                fh2.Position = [1100 80 1064 1064];
                ah2 = tight_subplot(subplotn(1),subplotn(2),[.01 .01],[0.01 0.01],[0.01 0.01]);

                stainmap = ones(256,3)*.05;
                stainmap(:,2) = linspace(.05,.9,256);
                unif_caxis = [prctile(vertcat(MPStats.stain_int),1) prctile(vertcat(MPStats.stain_int),99)];
            end

            if plotstain2
                if isfield(MPStats,[stainname2 '_int'])
                    % IO_Stain
                    fh3 = figure('Name',stainlabel);
                    fh3.Position = [1100 80 1064 1064];
                    ah3 = tight_subplot(subplotn(1),subplotn(2),[.01 .01],[0.01 0.01],[0.01 0.01]);

                    stainmap2 = ones(256,3)*.05;
                    stainmap2(:,1) = linspace(.05,.9,256);
                    stainmap2(:,3) = linspace(.05,.9,256);
                    unif_caxis2 = [prctile(vertcat(MPStats.([stainname2 '_int'])),1) prctile(vertcat(MPStats.([stainname2 '_int'])),99)];
                end
            end
        end
    end
    
    if exist('plotmask','var')
        if plotmask == 1           
            fh6 = figure('Name','Mask');
            fh6.Position = [1200 150 1064 1064];
            ah6 = tight_subplot(subplotn(1),subplotn(2),[.01 .01],[0.01 0.01],[0.01 0.01]);                       
        end
    end
    
    if exist('plotcurv','var')
        if plotcurv == 1
            
            fh4 = figure('Name','Mean Curvature deviation');
            fh4.Position = [1200 150 1064 1064];
            ah4 = tight_subplot(subplotn(1),subplotn(2),[.01 .01],[0.01 0.01],[0.01 0.01]);
            
            fh5 = figure('Name','Gaussian Curvature');
            fh5.Position = [1200 150 1064 1064];
            ah5 = tight_subplot(subplotn(1),subplotn(2),[.01 .01],[0.01 0.01],[0.01 0.01]);
            
        end
    end
           
    for iMP = (1+(jbatch-1)*ntot):min(ntot*jbatch,length(MPStats))
        [edgecoor,centroidcoor,edgecoor_sph,TRI_nosph] = check_processing_status(MPStats(iMP));
        ThetaPhiR = edgecoor_sph;
        R = ThetaPhiR(:,3);
        if mollweid_proj
            [Xmw,Ymw,TRI_nosph] = MollweidProj(ThetaPhiR(:,1),ThetaPhiR(:,2));
            ThetaPhiR(:,[1 2]) = [Xmw,Ymw];
        end
        
        % Plot the shape
        axes(ah(iMP-(jbatch-1)*ntot)); %#ok<LAXES>
        Make_2D_plot(TRI_nosph,ThetaPhiR,R);
        colormap(brewermap([],'RdBu'))
        caxis(median(R)+[-.6 .6])

        %Plot the LifeAct
        if exist('plotstain','var') 
            if plotstain
                if isfield(MPStats,[stainname '_int'])
                   axes(ah2(iMP-(jbatch-1)*ntot)); %#ok<LAXES>
                   stainint = MPStats(iMP).([stainname '_int']);
                   Make_2D_plot(TRI_nosph,[ThetaPhiR(:,1:2) stainint],stainint);
                   colormap(stainmap)
                   caxis([prctile(stainint,1) prctile(stainint,99)])
                   %caxis(unif_caxis)
                end

                if plotstain2
                    if isfield(MPStats,[stainname2 '_int'])
                        axes(ah3(iMP-(jbatch-1)*ntot)); %#ok<LAXES>
                        stainint = MPStats(iMP).([stainname2 '_int']);
                        Make_2D_plot(TRI_nosph,[ThetaPhiR(:,1:2) stainint],stainint);
                        colormap(stainmap2)
                        caxis([prctile(stainint,1) prctile(stainint,99)])
                        %caxis(unif_caxis2)
                    end
                end
            end
        end
        
        if exist('plotmask','var')
            if plotmask == 1
                %Plot the Mask
                set(groot,'CurrentFigure',fh6)
                fh6.CurrentAxes = ah6(iMP-(jbatch-1)*ntot); 
                stainint = MPStats(iMP).isincontact;
                Make_2D_plot(TRI_nosph,[ThetaPhiR(:,1:2) stainint],stainint);
                colormap(gray)
                caxis([0 1])
            end
        end     
        
        if exist('plotcurv','var')    
            if plotcurv == 1
                
                [ThetaPhiRsmooth] = Smooth_spherical_mesh(edgecoor_sph,sqrt(1^2/pi),'arc',1,'mean');
                [X,Y,Z] = sph2cart(ThetaPhiRsmooth(:,1),ThetaPhiRsmooth(:,2),ThetaPhiRsmooth(:,3));    
            
                FV.faces    = MPStats(iMP).TRI_Connectivity_SPHbound;
                FV.vertices = [X,Y,Z];
                pc = GetCurvatures(FV,0);
                gc = reshape(pc(1,:).*pc(2,:),[],1);
                mc = reshape(mean(pc),[],1);
                mc = mc - median(mc);

            axes(ah4(iMP-(jbatch-1)*ntot)); %#ok<LAXES>
            Make_2D_plot(TRI_nosph,[ThetaPhiR(:,1:2) mc],mc);
            colormap(greenmagenta)
            caxis([-0.3 0.3])

            axes(ah5(iMP-(jbatch-1)*ntot)); %#ok<LAXES>
            Make_2D_plot(TRI_nosph,[ThetaPhiR(:,1:2) gc],gc);
            colormap(turquoisebrown)
            caxis([-0.18 0.18])
                
            end
        end                              
    end
    
    AH{jbatch}  = ah ;
    if exist('plotstain','var') 
        if plotstain == 1
            AH2{jbatch} = ah2;
            if plotstain2
                if isfield(MPStats,[stainname2 '_int'])
                    AH3{jbatch} = ah3;
                end
            end
        end
    end
    if exist('plotcurv','var')    
        if plotcurv == 1
            AH4{jbatch} = ah4;
            AH5{jbatch} = ah5;
        end
    end
    if exist('plotmask','var')
        if plotmask == 1
            AH6{jbatch} = ah6;
        end
    end
    
    
end

if rem(iMP,ntot)~=0
    for iax = ntot:-1:rem(iMP,ntot)+1
        delete(ah(iax))
        ah(iax) = [];
        if exist('plotstain','var')
            if plotstain == 1
                if isfield(MPStats,[stainname '_int']) ;    delete(ah2(iax));   ah2(iax) = [];       end
                if plotstain2
                    if isfield(MPStats,[stainname2 '_int']);    delete(ah3(iax));   ah3(iax) = [];       end
                end
            end
        end
        if exist('plotmask','var')
            if plotmask == 1;   delete(ah6(iax));   ah6(iax) = [];      end
        end
        if exist('plotcurv','var')    
            if plotcurv == 1   
                delete(ah4(iax));   delete(ah5(iax));
                ah4(iax) = []   ;   ah5(iax) = []   ;
            end
        end
    end    
end

AH{jbatch}  = ah ;
if exist('plotstain','var') 
    if plotstain == 1
        AH2{jbatch} = ah2;
        if plotstain2
            if isfield(MPStats,[stainname2 '_int'])
                AH3{jbatch} = ah3;
            end
        end
    end
end
if exist('plotcurv','var')    
    if plotcurv == 1
        AH4{jbatch} = ah4;
        AH5{jbatch} = ah5;
    end
end
if exist('plotmask','var')
    if plotmask == 1
        AH6{jbatch} = ah6;
    end
end

%% Plot only the ROI-segmented mean curvature (for display purposes, 1-by-1)

% Open ROI structure file first
ThetaPhi_InROI = inpolygon(ThetaPhiR(:,1),ThetaPhiR(:,2),roi_position.position(:,1),roi_position.position(:,2));
mc_in = mc .* ThetaPhi_InROI;

% Plot actin map with ROI overlay
figure
Make_2D_plot(TRI_nosph,ThetaPhiR(:,[1 2]),stainint);
colormap(stainmap)
hold on, h = patch(roi_position.position(:,1),roi_position.position(:,2),[1 1 1],'FaceAlpha',0.2,'EdgeColor',[1 1 1],'LineWidth',2,'LineStyle','--');

figure
Make_2D_plot(TRI_nosph,ThetaPhiR(:,[1 2]),mc_in);
colormap(turquoisebrown) % sometimes I use this, sometimes greenmagenta
cmap = colormap;
cmap(129,:) = [1 1 1];
colormap(cmap)
caxis([-0.6 0.6])

%% Equalize axes and color scales

normalize_radius = 0;

if contains(get(get(ah(1),'parent'),'Name'),'Sphere')
    allxlim = vertcat(vertcat(AH{:}).XLim);
    allylim = vertcat(vertcat(AH{:}).YLim);
    allzlim = vertcat(vertcat(AH{:}).ZLim);
    lowerlims = min([allxlim(:,1), allylim(:,1), allzlim(:,1)]);
    upperlims = max([allxlim(:,2), allylim(:,2), allzlim(:,2)]);
    for i = 1:length(AH)
        [AH{i}.XLim] = deal([lowerlims(1) upperlims(1)]);
        [AH{i}.YLim] = deal([lowerlims(2) upperlims(2)]);
        [AH{i}.ZLim] = deal([lowerlims(3) upperlims(3)]);
    end
end

if normalize_radius
    for i = 1:length(AH)
        [AH{i}.CLim] = deal([ah(1).CLim]+[-.4 .4]);
    end
end

if exist('ah2','var')
    for i = 1:length(AH)
        [AH2{i}.CLim] = deal(prctile(vertcat(MPStats.stain_int),[1 99]));
        if contains(get(get(ah(1),'parent'),'Name'),'Sphere')
            [AH2{i}.ZLim] = deal([lowerlims(3) upperlims(3)]);
            [AH2{i}.XLim] = deal([lowerlims(1) upperlims(1)]);
            [AH2{i}.YLim] = deal([lowerlims(2) upperlims(2)]);
        end
    end
end

if exist('ah3','var')
    for i = 1:length(AH)
        %[AH3{i}.CLim] = deal(prctile(vertcat(MPStats.stain_int),[1 99]));
        if contains(get(get(ah(1),'parent'),'Name'),'Sphere')
            [AH3{i}.ZLim] = deal([lowerlims(3) upperlims(3)]);
            [AH3{i}.XLim] = deal([lowerlims(1) upperlims(1)]);
            [AH3{i}.YLim] = deal([lowerlims(2) upperlims(2)]);
        end
    end
end

if exist('ah4','var')
    for i = 1:length(AH)
        %[AH4{i}.CLim] = deal(prctile(vertcat(MPStats.stain_int),[1 99]));
        [AH4{i}.CLim] = deal([-.5 .5]);
        if contains(get(get(ah(1),'parent'),'Name'),'Sphere')
            [AH4{i}.ZLim] = deal([lowerlims(3) upperlims(3)]);
            [AH4{i}.XLim] = deal([lowerlims(1) upperlims(1)]);
            [AH4{i}.YLim] = deal([lowerlims(2) upperlims(2)]);
        end
    end
end

if exist('ah5','var')
    for i = 1:length(AH)
        [AH5{i}.CLim] = deal(prctile(vertcat(MPStats.stain_int),[1 99]));
        if contains(get(get(ah(1),'parent'),'Name'),'Sphere')
            [AH5{i}.ZLim] = deal([lowerlims(3) upperlims(3)]);
            [AH5{i}.XLim] = deal([lowerlims(1) upperlims(1)]);
            [AH5{i}.YLim] = deal([lowerlims(2) upperlims(2)]);
        end
    end
end

if exist('ah6','var')
    for i = 1:length(AH)
        [AH6{i}.CLim] = deal(prctile(vertcat(MPStats.stain_int),[1 99]));
        if contains(get(get(ah(1),'parent'),'Name'),'Sphere')
            [AH6{i}.ZLim] = deal([lowerlims(3) upperlims(3)]);
            [AH6{i}.XLim] = deal([lowerlims(1) upperlims(1)]);
            [AH6{i}.YLim] = deal([lowerlims(2) upperlims(2)]);
        end
    end
end

%% Save images of MC plots one by one for stacking

for iMP=1:length(MPStats)
    figure
            [edgecoor,centroidcoor,edgecoor_sph,TRI_nosph] = check_processing_status(MPStats(iMP));
            % Calculate curvatures.
                theta = edgecoor_sph(:,1); phi = edgecoor_sph(:,2);
                ThetaPhiR = edgecoor_sph;         
                [ThetaPhiRsmooth] = Smooth_spherical_mesh(edgecoor_sph,sqrt(1^2/pi),'arc',1,'mean');
                [X,Y,Z] = sph2cart(ThetaPhiRsmooth(:,1),ThetaPhiRsmooth(:,2),ThetaPhiRsmooth(:,3));    
                FV.faces    = MPStats(iMP).TRI_Connectivity_SPHbound;
                FV.vertices = [X,Y,Z];
                pc = GetCurvatures(FV,0);
                gc = pc(1,:).*pc(2,:);
                mc = mean(pc);
            % Plot curvature map
                ph_synapse = patch('faces',TRI_nosph,'vertices',[ThetaPhiR(:,1:2) mc'],'facevertexcdata',mc','facecolor','flat','edgecolor','none');
                colormap(turquoisebrown)
                caxis([-1 1])
                % Alternative scaling: caxis([prctile(mc,1) prctile(mc,99)]) 
                box off
                axis off
end

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


% Set 2D plot options
function[ph] = Make_2D_plot(TRI,vert,color)

    ph = patch('faces',TRI,'vertices',vert,'facevertexcdata',mean(color(TRI),2),'facecolor','flat','edgecolor','none');
    view(3)
    axis off
    axis equal
    view(0,90);
    set(gca,'Clipping','off')

end

