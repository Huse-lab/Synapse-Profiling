% ImageAnalysis_SingleParticleMovie

% For a general overview of the methodology, see:
% Vorselen, D. et al. Microparticle traction force microscopy reveals subcellular force exertion patterns in immune cell?target interactions. Nat. Commun. 11, 20 (2020).
% (https://www.nature.com/articles/s41467-019-13804-z)

% This code will allow you to analyze the deformation of a single particle
% over time or multiple separate particles, as long as each image only contains one particle. 

% Run one section at a time by pressing ctrl-enter (or cmnd-enter) 

% If only a fluorescent signal from a particle (i.e. no cell-stain or label 
% for the free particle surface), the last useful cell is "Detect edges using 
% the 3D sobel operator and superlocalize edges using gaussian fitting". 
% MPStats (the structure containing the data) and text files with the edge 
% coordinates and triangulation, can still be saved using the last two cells

% Useful_plots.m can be used to make a couple of example plots

% If you use any part of this code , please cite the manuscript above in any publications

%% Read the data files

% Initialize a parallel pool and add the required dynamic javapaths
p = gcp;
w = warning('off', 'MATLAB:Java:DuplicateClass');
javaaddpath(fileparts(which('ParforProgressMonitor.class')));
warning(w);
clear('p','w');

% Open UI to select the data files
[data,~,FileName,PathName] = ReadBFImages('.tif');

%% Extract the images and required metadata

if exist('OGdata','var'); data = OGdata;  end
[MPStats,data1                ] = Extract_Essential_Metadata(data,FileName,PathName,'istimeseries',1);
[IM3Dcell,LLSMMask,zcorrfactor] = Extract_Images(data1,MPStats,MPStats(1).MPchannel);
[MPStats.IM3D                 ] = IM3Dcell{:};
[MPStats.LLSMMask             ] = LLSMMask{:};
[MPStats.PixelSizeZ_original  ] = MPStats.PixelSizeZ;
[MPStats.PixelSizeZ           ] = ArraytoCSL([MPStats.PixelSizeZ]/zcorrfactor);
[MPStats.zcorrFactor          ] = deal(zcorrfactor);
  
MPStats = orderfields(MPStats);
if ~all(size(data)==size(data1));     OGdata = data;     data = data1;       end        
clear('IM3Dcell','LLSMMask','zcorrfactor','data1');

%% Threshold the images and identify particles

% Increasing the "UseIncreasedThreshold" value (between 0 and 1) excludes less particles that are
% close to/touching the border. Warning: this could result in errors later on. Using the watershed
% option increases computational time by a lot, but allows analyzing movies with adjacent particles.
Opts = {'UseIncreasedThreshold',0,'watershed',1};

MPStats = Threshold_images_and_identify_particles(MPStats,Opts);
clear('Opts');

%% Superlocalize particle edges and triangulate surface

% Ask the user about the use of sobel filtering and how much smoothing to apply to the image
edgedetectionsettings = inputdlg({'Use Sobel filter ("sobel") or derivative of line profiles ("direct")?',...
    'How much smoothing do you want to apply (in pixels)?',['Do you want to use stable mode? '...
    '(increases computational complexity and time, but can work better especially with adjacent particles)'],...
    'Desired (max) spacing between the points (in nm)?',['Do you want to use a regular spaced grid (reg) or '...
    'equidistant points (equi)?']},'Edge Detection Settings',[1 80],{'direct','1','0','250','equi'});

% Apply sobel filtering (if desired)
usesobel = strcmpi(edgedetectionsettings{1},'sobel');
if usesobel
    XYZedges = Apply_3DSobel_filter({MPStats.IM3D},edgedetectionsettings);
    [MPStats.IMedges] = XYZedges{:};
end

[MPStats,Residuals_Problems_excluded] = Superlocalize_edges_using_gaussian_fitting(MPStats,...
    str2double(edgedetectionsettings{4})/1000,edgedetectionsettings{5},'DirectDerivative',~usesobel,...
    'usestablemode',str2double(edgedetectionsettings{3}),'smooth',str2double(edgedetectionsettings{2}),'RadialBounds',[0.25 0.7]);

% Triangulate surface and determine particle statistics
MPStats = Triangulate_surface_and_determine_particle_statistics(MPStats);

clear('XYZedges','edgedetectionsettings','usesobel')

%% Determine particle coverage by a secondary signal

if exist('data','var')
    Stain_channel = Find_Stain_Channel(data(1),MPStats);
    MPStats = Analyze_Secondary_Signal(data(:,1),MPStats,[],Stain_channel,'stain');
else
    MPStats = Analyze_Secondary_Signal([],MPStats,[],[],'stain');
end

clear('Stain_channel')

%% Convert secondary signal to mask and align cups

% Choose to use the maximum intensity or integrated intensity. To choose, you can plot both by running this line:
% figure; subplot(2,1,1); imagesc(max(MPStats(1).IMstain_radial,[],3)'); axis equal; axis off; subplot(2,1,2); imagesc(sum(MPStats(1).IMstain_radial,3)'); axis equal; axis off;
% The upper one is max intensity, lower integrated intensity
userinput = inputdlg('Use maximum intensity (max) or integrated (int) intensity signal?',...
            'Stain analysis options',[1 70],{'int'});
use_integrated_intensity = strcmp(userinput,'int');
stain_indicates_contact  = 1;

MPStats = Convert_secondary_signal_to_mask(MPStats,stain_indicates_contact,'use_integrated_intensity',use_integrated_intensity,'use_global_threshold',0);
MPStats = Determine_base_position_and_align(MPStats,'BaseLat',-pi/2,'BaseColongitude',0);

clear('use_integrated_intensity','stain_indicates_contact')

%% Optional: determine particle coverage by additional signals

% Check if at least one stain has been analyzed so far
if ~isfield(MPStats,'IMstain')
    error(['If this is the first fluorescent signal you analyze, use "Determine'...
        ' particle coverage by a secondary signal" two cells up'])
end

if exist('data','var')
    [MPStats,IMStats] = Analyze_Secondary_Signal(data,MPStats,IMStats,'stain',1);
else
    [MPStats,IMStats] = Analyze_Secondary_Signal([],MPStats,IMStats,'stain',1);
end

%% Copy_Change_Mask_MDJ: Took apart Daan's old version of the code (pre-2021-07-08) to be compatible with new version.

warning('off','all');

for idx=1:length(MPStats)
    if idx==1
        MPStats_dummy = Update_Mask_GUI_MDJ(MPStats(idx));
        MPStats_dummy.FileName
    elseif idx>1
        MPStats_dummy(idx) = Update_Mask_GUI_MDJ(MPStats(idx));
        MPStats_dummy(idx).FileName
    end
end

clear MPStats;
MPStats = MPStats_dummy;
% clear MPStats_dummy

disp('Done updating all masks for this MPStats file.');
