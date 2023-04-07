# Synapse-Profiling
Analysis Pipeline accompany manuscript: de Jesus &amp; Settle et al 2023

DAAM Particle Image Processing and Synapse Profiling Guide

This repository contains functions, scripts, and test-cases to demonstrate and reproduce the topographical analysis of immune synapses interacting with DAAM particles, as described in the associated manuscript: 

"Topographical analysis of immune cell interactions reveals a biomechanical signature for immune cytolysis"
Miguel de Jesus*, Alexander Settle*, Daan Vorselen, Michael Galiano, Yung Yu Wong, Tian-Ming Fu, Endi Santosa, Fella Tamzalit, Mitchell S. Wang, Zhirong Bao, Joseph C. Sun, Pavak Shah, Julie Theriot, and Morgan Huse†.
In Preparation 2023 [Title, Author list, Citation subject to change after submission] 

* equal contribution
† corresponding author


Contents Overview
1. Preprocessing Pipeline: Convert raw image files into MPStats, annotate ROIs, and convert to MPRender for faster processing
2. Synapse Analysis Pipeline: extract useful information from MPRender and ROI and store as a synapseResults MATLAB object for batching and comparing across particles
3. Example .tifs, MPRenders, and ROI files that can be run through the pipeline. 
4. Functions for all computations in the pipeline
5. Dependencies: non-original functions and useful code from other sources 


Requirements: MATLAB 2021 or later.
Set MATLAB path to contain all subfolders within this repository Synapse-Profiling-HuseLab. All dependencies and functions should be contained within. Report any missing links or bugs to alex.h.settle@gmail.com or husem@mskcc.org


PREPROCESSING
FROM IMAGE TO MPRENDER & ROI FILES
The pre-processing step is mostly adapted from Vorselen, D. et al. Microparticle traction force microscopy reveals subcellular force exertion patterns in immune cell target interactions. Nat. Commun. 11, 20 (2020). (https://www.nature.com/articles/s41467-019-13804-z)
See original paper for a general overview of the methodology. 

Step 1: Load .tiff stacks into a working folder
	- data should be two-color Z-stacks saved as one file per time-point as a .tif. Current version only works on .tifs 

Step 2: Convert tif data into MPStats. 
		- Open ImageAnalysis_tiftoMPStats.m, run each section individually:
	•	%% Read the data files
	•	%% Extract the images and required metadata
	•	Click on blob channel
	•	Z-correction factor "1"
	•	If needed, enter pixel size information: for Vt-iSIM, 0.108 x 0.108 x 0.300 nm^3
	•	%% Threshold the images and identify particles
	•	Particles that encounter problems at this step are excluded from analysis. 
	•	%% Superlocalize particle edges and triangulate surface
	•	Pick smoothing of "0" or "1"; for most analyses in the paper 0 was used (Fig. 1 demonstrations use 1 for clarity)
	◦	%% Determine particle coverage by a secondary signal
	•	No changes from default parameters
	•	%% Convert secondary signal tomask and align cups
	•	Stain analysis options: "max"
	•	%% Realign cup
	•	Realign the particle such that actin mask is on the center-left third of the theta-phi map
	•	If masking is done properly, simply clicking "Realign" achieves this; otherwise, click "Clear Mask(!)" re-draw the mask manually, and realign.
	•	Clicking "Done" updates MPStats structure
	•	%% Save MPStats file
	•	Assign a name to the output file. This will contain all the selected data.
	•	For manual inspection of particle renderings: run generatePlotsBatch.m
	•	Saves figures as .figs and .png 
	•	
Step 3: Reduce data size to save memory and processing time
	•	Open convert_MPStatstoMPRender.m
	•	run Select Folder of Interest, which should contain any number of MPStats files
	•	run Loop section, this will automatically run through everything and convert to MPRender
	•	MPRender contains all of the information in MPStats except for raw image data. 
	•	This step is technically optional but it will save a lot of time and computational power so it is highly recommended.


Step 4: Manually Select ROI by drawing on 2D projections. 
	⁃	To allow consistent comparisons between particles with different secondary signals that may threshold differently, we chose to manually annotate regions of interest (ROI) for all particles. An automated version of this is in progress, but for de Jesus and Settle et al, all ROIs were manually drawn by enclosing the contiguous area of secondary signal. 
	⁃	Open and run selectROI_batch.m
	⁃	if cell stain exists, chooses option 1, otherwise 0
	⁃	This script will loop through all the MPStats in the folder and allow you to draw each one at a time.




SYNAPSE ANALYSIS PIPELINE
In this section, we take the MPStats particle data and the ROI information and calculate all relevant statistics calculated in the manuscript (deJesus, Settle et al. 2023) and store them in a large MATLAB struct called synapseResults. This struct can either be a single 1x1 struct or a struct array, depending on the size of the MPStats and ROI inputs. See below for detailed description of the synapseResult fields.

To generate results for a single synapse:
	- Run ProfileSingleSynapse.m
		- The script will prompt you to select the MPStats file and corresponding ROI file via a GUI.
		- The script will save the synapseResults struct as a .mat file and a .png showing a overview of relevant statistics calculated and their associated data

To generate results for a large batch of data
	- Transfer all MPRenders and ROIs into one directory
		- Naming conventions are important, make sure the ROI file names match the MPRenders (the preprocessing will output them according to the convention)
	- Run ProfileBatchSynapses.m
		- select the directory of interest
		- this will take ~10-30 seconds per synapse depending on your computer
		- The script will print each filename as it calculates so you can track progress, set verbose = 0 to turn this off. 
		- MPStat pools will be kept together as one synapseResults struct, but separate ones will be separated. 
	- This will also ask whether you want to make a overview/results summary page for each synapse. This adds time/memory to the total process, and is unnecessary for downstream processing. You can also generate these later using generateSynapseResultsPage.m 


Synapse Results Data Structure
- FileName: string, the name of the original .tif file analyzed
- MeanCurvatureGrid: 50x50 double, a gridded matrix of the synapse ROI, with values at each position defined by the mean curvature at that location
- GaussianCurvatureGrid: 50x50 double, a gridded matrix of the synapse ROI, with values at each position defined by the gaussian curvature at that location
- StainIntensityGrid: 50x50 double, a gridded matrix of the synapse ROI, with values at each position defined by intensity of fluorescence of the cell stain / label
- gridScale_um: scalar, the width/height of the ROI in µm in cartesian space. Because the grids are projected onto a normalized grid, this factor can be used to back-calculate feature sizes in real space
- ZernikeModalCoefficients: 136x1 double, modal coefficients for each Zernike Polynomial computed up to 15 radial degrees, with all azimuthal orders included. The array is ordered by the OSA/ANSI standard indices convention. See Materials and Methods in the associated manuscript. 
- ZernikeRDCoefficients: 1x72 double, Rotationally degenerate coefficients transformed from the modal coefficients. See Methods for details on this methodology.
- ZernikePatternComplexity: integer, the Pattern Complexity Number, as defined by the degrees of Zernike Polynomial computations required to faithfully reconstruct the original data to within 0.8 Pearson correlation. 
- ZernikeComplexityTrace: 1x30 double, array detailing the Pearson correlation between the original data and reconstructed data with increasing Zernike computations, from 1 radial degrees (1 function) to 30 radial degrees. Mostly useful for plotting.
- RadialCurvatureData: MATLAB struct containing Nx3 arrays defining the edge coordinates within the ROI, as well as their associated curvatures and their normalized distance to the ROI center. These values are used for calculating radial curvature profiles.
- MeanCurvatureRadialProfile: 20 x 10 double, binned RadialCurvatureData for mean curvature, each column is a histogram of mean curvature values within a radial bin. The column binning spans from 0->1 with 0.1 iintervals, and the row binding is -1->1 with 0.1 intervals. See Materials and Methods for more on the methodology.
- GaussianCurvatureRadialProfile: 20 x 10 double, binned RadialCurvatureData for gaussian curvature, each column is a histogram of mean curvature values within a radial bin. The column binning spans from 0->1 with 0.1 iintervals, and the row binning is -0.3->0.3 with 0.03 intervals. See Materials and Methods for more on the methodology
- BulkDeformationResults: MATLAB struct containing data on the bulk 3D-deformation of the particle. See function calculateBulkDeformation documentation for complete list of calculated items.
- ProtrusionStats: MATLAB struct containinig data on all negative-curvature features, see featurePicker for details.
- PeakStats: MATLAB struct containinig data on all positive-curvature features, see featurePicker for details.
- FeatureSummary: MATLAB struct containing summarized/averaged information on protrusions and peaks for a given synapse.




