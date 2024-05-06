% Zernike Complexity as RECONSTRUCTION RESISTANCE
% Uploaded 12/9/21
% Modified by MDJ 2022-06-15

% Modification by MDJ 20231129: Try adopting as alt. pattern complexity #
% the integrated Pearson fits

% Modification by MDJ 20240108: Discussed a new method to calculate
% Z-complexity with AHS, as an area instead of a threshold cut-off. We
% measure complexity as the area OVER the curve, bounded by a line at
% y=R=1. This means that the simpler patterns reach this area much quicker
% than complex ones. So we turn it into a "forward" complexity metric by
% subtracting it from the maximum area, which is 1 x max_order = max_order.
% i.e., A rectangle with 1 as its maximum height and max_order as its maximum length.

%% Import grids into a cell array 
% This code segment will search a folder containing grids and import
% them into 50x50 data structures and into an array. Note that this
% will ask you to import both MC and GC regardless of whether you want to
% include both or either in your analysis. This is just to simplify the
% code and not go crazy.

fprintf('Please select folder of interest containing MC grids \n')
mc_file_path = uigetdir(pwd,'Please select folder of interest containing MC grids');
% fprintf('Please select folder of interest containing GC grids \n')
% gc_file_path = uigetdir('Select Folder of interest');

%Set file paths based on user input
current_dir = pwd;
addpath(mc_file_path);
files = dir(mc_file_path);

%Initialize
% gc_synapses = {};
mc_synapses = {};
names = {};
labels = {};
colors = [];

for file = files'
    if contains(file.name,'_MeanC') % check for only 50x50 grid files
%         load(file.name)  %loaded synapse will have variable name 'mc_i'
%         cd(gc_file_path)
        load(file.name)
%         addpath(mc_file_path)
        
%         gc_synapses = [gc_synapses; gc_i];
        mc_synapses = [mc_synapses; mc_i];
        names = [names;file.name]; %store file names for plotting/organization
                
    end
end

%% Set Parameters

max_order = 40;

%% Calculate Pattern Complexity numbers 

L = size(mc_i,1);
X = -1:2/(L-1):1;
[x,y] = meshgrid(X);
x = x(:); y = y(:);
[theta,r] = cart2pol(x,y);

pearson_array = struct('R',[]);

for i = 1:length(mc_synapses)
    disp(i);
    mc_i = mc_synapses{i};
    fits = [];
    reconstructions = {};

    for j = 1:max_order
        N = []; M = [];
        for n = 0:j
            N = [N n*ones(1,n+1)];
            M = [M -n:2:n];
        end
        is_in_circle = ( r <= 1 );
        Z = zernfun(N,M,r(is_in_circle),theta(is_in_circle));
        a = Z\mc_i(is_in_circle);
        mc_i(~is_in_circle) = 0;
        reconstructed = zeros(size(mc_i));
        reconstructed(is_in_circle) = Z*a;
        reconstructions = [reconstructions; reconstructed];
        fit = corr2(mc_i,reconstructed);
        fits = [fits; fit];
    end
    
    pearson_array(i).R = fits;

end

%% Restructure before plotting

complex_struct = struct('FileName',{},'Complexity_resistance',[]);

for k = 1:length(names)
    complex_struct(k).FileName = names(k);
    if exist('pearson_array')
        area_under_curve = trapz(pearson_array(k).R);
        complex_struct(k).Complexity_resistance = max_order - area_under_curve;
    end
end
complex_struct = assign_celltype(complex_struct);
complex_struct = sortStruct(complex_struct,'Celltype');

%% Perform all  
adjustment of Celltypes AND RESORTING before proceeding

%% Assign labels
cell_labels = [{complex_struct.Celltype}];
labels = unique(cell_labels);

%% Plot Histogram

figure
hold on
for k = 1:length(labels)
    histogram([complex_struct(strcmp({complex_struct(:).Celltype},labels{k})).Complexity_resistance],1:35,'Normalization','pdf','DisplayStyle','stairs')
    xlabel('Zernike Complexity')
    ylabel('pdf')
end
legend(labels)

% saveas(gcf,'IMG117_complexity_cutoff75.png')

%% Correlation plot of PC1 vs. Settle Complexity

figure
% Conduct PCA analysis on a pre-set coeff_struct.mat file, with the EXACT
% matched celltype labels & order (Sorted after manually annotated).
    feature_vect = [];
    celltype_vect = {};
    % Fetch data
    for i = 1:length(coeff_struct_mc)
        celltype_vect = [celltype_vect; coeff_struct_mc(i).Celltype];
%       if strcmp(Grid_type,'mc')
            feature_vect = [feature_vect; coeff_struct_mc(i).C_adjusted_zernike];
%        elseif strcmp(Grid_type,'gc')
%            feature_vect = [feature_vect; complex_struct(i).GC_adjusted_zernike];
%        elseif strcmp(Grid_type,'both')
%            feature_vect = [feature_vect; [complex_struct(i).MC_adjusted_zernike complex_struct(i).GC_adjusted_zernike]];
%        end
    end
[coeff, score, latent, tsquared, explained] = pca(feature_vect);

scatter(score(:,1),[complex_struct.Complexity_resistance],15,[0 0.6 0.8],'filled')
pbaspect([1 1 1])
xlabel(['Principal Component 1 (',num2str(explained(1)),' %)'])
ylabel('Synapse complexity')

x = [score(:,1)];
y = [complex_struct.Complexity_resistance]';
% b = x\y;
hold on
yfit = x*b;
title(['Zernike complexity analysis (n = ' num2str(length(cell_labels)) ')'])
plot(score(:,1),yfit,'-','LineWidth',3,'Color',[.85 .25 .25])

R_fit = corr2(x,y);
v = axis;
text(v(2)*0.6,v(4)*0.95,['R = ' num2str(R_fit)])

%% Correlate Complexity (complex_struct) vs. dV (LPA_Stats3D)

if ~exist('complex_struct','var')
    error('Complexity structure is not loaded into the file.')
elseif ~exist('LPA_Stats3D','var')
    error('Strength data is not loaded into the file.')
end

% Check even sorting of files

matchcount = [];
% complex_struct = sortStruct(complex_struct,'Celltype');
% LPA_Stats = sortStruct(LPA_Stats,'Celltype');

names_complex = {complex_struct.FileName};
names_indentation = {LPA_Stats3D.FileName};
for idx=1:length(names_complex)
    trunc = extractBefore(names_indentation{idx},'.');   
    matchcheck = or(contains(char(names_complex{idx}),trunc),contains(trunc,char(names_complex{idx})));
    matchcount = [matchcount matchcheck];
end

if any(matchcount == 0)
    % Re-sort both by filename:
    dummy_complex = struct2table(complex_struct);
    sorted = sortrows(dummy_complex,'FileName');
    complex_struct = table2struct(sorted);
    dummy_complex = struct2table(LPA_Stats3D);
    sorted = sortrows(dummy_complex,'FileName');
    LPA_Stats3D = table2struct(sorted);    
end

% There is another way as well to do this: 
% idx_mismatch = ismember(extractBefore(names_complex,'_50x50_'),extractBefore(names_indentation,'.tif'));
% Then delete the mismatch == 0. Haven't coded this in yet.

clear matchcount names_complex names_indentation trunc matchcheck dummy_complex sorted

% Plot

x_strength = [LPA_Stats3D.Volume_Deviation_ROI];
y_complexity = [complex_struct.Complexity_resistance];
ref_celltype = {complex_struct.Celltype};
[list_celltype] = unique(ref_celltype);

% figure
    x_bounds = [min(x_strength) max(x_strength)];
    y_bounds = [min(y_complexity) max(y_complexity)];
for idx=1:length(list_celltype)
    % Find indices per data category (celltype)
%     n_col = ceil(sqrt(length(list_celltype)));
%     n_row = ceil(length(list_celltype)/n_col);
    figure
    [within_celltype] = ismember(ref_celltype,list_celltype{idx});
    % Plot scatter gallery
%     subplot(n_row,n_col,idx), hold on
    scatter(x_strength(within_celltype),y_complexity(within_celltype),15,kulay(idx,:),'o',"filled");
    xlabel('Indentation strength');ylabel('Spectral complexity');    
%     xlim(x_bounds);ylim(y_bounds); 
    xlim(x_bounds); ylim(y_bounds); pbaspect([1.5 1 1]);
     title(list_celltype{idx},'Interpreter','none');
    R = corr2(x_strength(within_celltype),y_complexity(within_celltype));
    text(2,y_bounds(2)+0.25,strcat("r = ",num2str(R)));
end
