% MPStats_to_MPRender_v6.m
%
% To enable higher-throughput and faster analysis of the data, this script
% reduces the MPStats files by removing raw data and leaving only the
% rendered data, reducing the file size significantly, creating a new file
% labeled as MPRender. MPRender will be saved with 'MPStats' as the
% variable name to stay consistent and compatible with the rest of my code
% Alex Settle & Miguel de Jesus
% Memorial Sloan Kettering Cancer Center
% Morgan Huse Laboratory, Department of Immunology
% 2023


% 2020-12-14, AHS: to render stain data in the edge coordinate space, thus
% allowing IMstain to be deleted. 

% 2021-03-30, MDJ: Added -v7.3 .mat file saving at the end.

% 2022-05-25, MDJ: By testing it without any of the Determine_Coverage
% steps, I figured the new Daanalysis code could run fine and go through
% with the correct stain coordinate system. So deleted the
% Determine_Coverage steps.

%% Select Folder of interest
fprintf('Please select folder of interest containing MPStats \n')
file_path = uigetdir('Select Folder of interest');
file_list = dir(file_path);

%% Loop through each MPStats file in the folder, remove raw data and save MPRender

for i = 1:length(file_list)
   if isempty(strfind(file_list(i).name,'MPStats_'))
       continue
   end
   disp(file_list(i).name)  
   file_name = fullfile(file_path, file_list(i).name);
   render_filename = strrep(file_name,'MPStats','MPRender');
   
   %check if render file already exists
   if isfile(render_filename)
       overwrite = questdlg('A render.mat file already exists for this image. Would you like to overwrite?', ... 
           file_list(i).name,'Yes','No','No');
       if strcmp(overwrite,'No')
           continue
       end
   end
   
   load(file_name);
   
   %Convert stain coverage data to match coordinates - THIS PART HAS BEEN
   %DELETED. REFER TO _v3 from the OBSOLETE folder for the reference code.

   if isfield(MPStats,'IM3D')
       MPStats = rmfield(MPStats,'IM3D');
   end
   if isfield(MPStats,'IMedges')
       MPStats = rmfield(MPStats,'IMedges');
   end
   if isfield(MPStats, 'IMlm')
       MPStats = rmfield(MPStats,'IMlm');
   end
   if isfield(MPStats, 'IM_bgcorr')
       MPStats = rmfield(MPStats,'IM_bgcorr');
   end
   if isfield(MPStats, 'IMstain')
       MPStats = rmfield(MPStats,'IMstain');
   end
   if isfield(MPStats, 'IMstain_bgcorr')
       MPStats = rmfield(MPStats,'IMstain_bgcorr');
   end
   if isfield(MPStats, 'IMstain_radial')
       MPStats = rmfield(MPStats,'IMstain_radial');
   end
   
   if ~contains(render_filename,'MPRender') % doublecheck file
       printf('Error, trying to overwrite wrong file')
   else
       fprintf(strcat('Saving','->',render_filename,'\n'))
       save(render_filename,'MPStats','-v7.3')
   end
   
   close
   
end

disp('Done converting MPStats to MPRenders.')
 