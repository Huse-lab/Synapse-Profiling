%% Autocorrelation through time (Kinetic analysis of grids)
% 2023-01-18: MDJ: Upon Alex's advice for Morgan's request on early
% supplementary figures

%% Initialize structure from grid data

file_path = uigetdir(pwd,"Open the MC grids folder.");
file_list = dir(fullfile(file_path,'*.mat'));
addpath(file_path);

grids = struct('filename',{},'mc_griddata',[]);
grids = timestamper(grids,'filename');

%% Auto-correlations of absolute topographies from grid data

    R_autocorr_topo = kinetic_autocorr(grids,'mc_griddata');
    figure, plot(R_autocorr_topo);
    
%% Auto-correlations of spatial freq. spectra from coeff_struct Zernike structure

if exist('coeff_struct','var')
    coeff_struct = timestamper(coeff_struct,'FileName');
end

% Perform auto-correlations on both MC_zernike & MC_adjusted_zernike

    R_autocorr_MCzernike = kinetic_autocorr(coeff_struct,'MC_zernike');
    R_autocorr_MCadjzernike = kinetic_autocorr(coeff_struct,'MC_adjusted_zernike');

    figure, hold on
    plot(R_autocorr_MCzernike), plot(R_autocorr_MCadjzernike);
    
%% Auto-correlations of pattern complexity

if exist('complex_struct','var')
    complex_struct = timestamper(complex_struct,'FileName');
end

% Perform auto-correlations on both MC_zernike & MC_adjusted_zernike

    R_autocorr_Complexity = kinetic_autocorr(complex_struct,'Complexity');

    figure, hold on
    plot(R_autocorr_Complexity);
    
%% Functions

function output_struct = timestamper(input_struct,fieldname)

        for idx=1:length(input_struct)
            timestamp = extractBetween(input_struct(idx).(fieldname),'_t','_50');
            timestamp = regexp(timestamp,'\d*','match');
            timestamp = str2double(timestamp{1});
            input_struct(idx).timestamp = timestamp; clear timestamp;
        end
        output_struct = sortStruct(input_struct,'timestamp');

end

function autocorr_data = kinetic_autocorr(input_struct,fieldname)

    % Create origial cell array to store the vector data.
    vect_orig = {input_struct.(fieldname)};

    % Prepare reference index list.
    vect_idx = 1:length(input_struct);
    vect_ref = cell2mat(vect_orig(vect_idx));

    autocorr_data = nan([1 size(input_struct,1)+1]);
    count = 0;
    % Use a rolling index array to create a dataset-wide query vector, then
    % perform dataset-wide autocorrelation between query & ref vectors.
    for idx=ceil(-0.5*length(vect_idx)):ceil(0.5*length(vect_idx))
        count = count+1;
        shuffle_idx = circshift(vect_idx,idx);
        test = vect_orig(shuffle_idx);
        vect_query = cell2mat(test);
        R_auto = corr2(vect_ref,vect_query);
        autocorr_data(count) = R_auto;
    end

end

