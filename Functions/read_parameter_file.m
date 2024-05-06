function parsed_data = read_parameter_file(param_file)

parsed_data.FileName = param_file;

meta_array = readmatrix(param_file);

%get number of protrusions
parsed_data.NumProtrusions = meta_array(1);

%get cluster size
parsed_data.MeanClusterDiameter = meta_array(2);
parsed_data.ClusterStdDev = meta_array(3);

%get bending rigidity
parsed_data.BendingRigidity = meta_array(4);


%get protrusion stiffness
parsed_data.ProtrusionStiffness = meta_array(5);

%get surface tension
parsed_data.SurfaceTension = meta_array(6);



end 
