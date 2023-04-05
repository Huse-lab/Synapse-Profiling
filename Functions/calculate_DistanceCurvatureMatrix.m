function mat_trans = calculate_DistanceCurvatureMatrix(data_input,nbin,cbin,curv_type)

 % Generate linearized distance and curvature arrays from the
 % cell-type-specific data structures.
        distances = data_input.Distance_NormalizedtoEdge(:)';
        if strcmp(curv_type,'mc')
            curvatures = data_input.MeanCurvature(:)';
        elseif strcmp(curv_type,'gc')
            curvatures = data_input.GaussianCurvature(:)';
        else
            error("please specify curvature type as 'mc' or 'gc'")
        end
    % Bin the arrays by radial distance. Define max as 1.0 to normalize
    % across cells.
        edges = linspace(0,1,(nbin+1));
        curvatures(2,:) = discretize(distances,edges);
    % Bin by curvatures
        curvature_window = struct('Radial_bin',[],'Curvatures',[],'Curvature_bin_freqs',[]);
        for idx=1:nbin
            curvature_window(idx).Radial_bin = idx;
                c_inShell = curvatures(1,curvatures(2,:)==idx);
                % Calculate frequencies of curvatures per bin/shell
                c_inShell(2,:) = discretize(c_inShell,cbin);
                c_binned = [1:length(cbin)-1];
                c_freq = histc(c_inShell(2,:),c_binned);
            curvature_window(idx).Curvatures = c_inShell(1,:);
            curvature_window(idx).Curvature_bin_freqs = [c_freq];
        end
        %clear c_inShell c_binned c_freq
    % Distance-curvature matrix (Location-strength binning)
        c_length = length(curvature_window(1).Curvature_bin_freqs);
        mat=zeros(nbin,c_length);
            for idx = 1:nbin
                for jdx = 1:c_length
                    mat(idx,jdx) = curvature_window(idx).Curvature_bin_freqs(jdx);
                end
            end
        mat_norm=zeros(nbin,c_length);
            for idx=1:nbin
                mat_norm(idx,:)=mat(idx,:)/sum(mat(idx,:));
            end
        mat_norm(isnan(mat_norm(:,:)))=0;
        mat_trans=mat_norm';
        
end