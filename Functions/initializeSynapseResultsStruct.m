function blankStruct = initializeSynapseResultsStruct()

blankStruct = struct('FileName',{},...
                 'MeanCurvatureGrid',[],...
             'GaussianCurvatureGrid',[],...
                'StainIntensityGrid',[],...
                      'gridScale_um',[],...
          'ZernikeModalCoefficients',[],...
             'ZernikeRDCoefficients',[],...
          'ZernikePatternComplexity',[],...
            'ZernikeComplexityTrace',[],...
               'RadialCurvatureData',{},...
        'MeanCurvatureRadialProfile',[],...
    'GaussianCurvatureRadialProfile',[],...
            'BulkDeformationResults',{});


end