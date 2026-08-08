# Gruzdeva-et-al.
This repository contains the code used to produce the analyses and figures in Gruzdeva et al., Hunger neurons track available food locations during foraging and spatial memory recall (2026).

The code processes fiber photometry recordings of AgRP neuron activity acquired alongside video tracking during three behavioral assays, and reproduces the quantifications, statistics, and figure panels reported in the manuscript.

**Requirements**
MATLAB R2019b or later
Statistics and Machine Learning Toolbox (fitlm, anova1, anovan, fitrm, ranova, multcompare, ranksum, kruskalwallis, ksdensity)
Image Processing Toolbox (imgaussfilt, used for spatial smoothing of occupancy and activity heatmaps)

No installation is required. Clone or download the repository, then add common_functions together with the module you intend to run, including its functions subfolder, to the MATLAB path:

addpath(genpath('common_functions'));
addpath(genpath('one_source_maze'));   % or 'three_sources_maze', or 'cheeseboard' 

**Repository organization**

The repository is divided into three modules, one per behavioral assay. Each contains a preprocessing script, a top-level pipeline script, and a functions folder holding the analysis and plotting routines it calls. Helpers shared across modules live in common_functions.

**Folders**

One_source_maze: 	Three-arm maze with a single food source; includes the GLM analyses
Three_sources_maze: 	Three-arm maze with three sources behind mesh doors, across sessions and light/dark conditions
Cheeseboard: 	Open cheeseboard arena, spatial learning, and memory recall
Common_functions: 	Routines used by more than one module


common_functions contains detect_food_runs, which segments continuous tracking into individual approach runs toward and away from a food source; findnearest, which returns the index of the sample closest to a target time and is used throughout preprocessing to align event times to tracking frames; and blueWhiteRed, the diverging colormap used for the activity and difference heatmaps. detect_food_runs applies the same detection logic to the single-source and three-source datasets, so run-based measures are directly comparable across the two paradigms; it accepts either data layout and is configured through an optional cfg struct documented in its header.

Running the code

Each module follows the same two-stage structure.

1. Preprocessing. The preprocessing scripts synchronize the raw acquisition streams into a single per-session matrix: photometry and analog TTL traces exported from the Doric console, centroid and DeepLabCut body-part coordinates from the top-view camera, and food-camera frame times. They compute dF/F from the 465 nm and 405 nm channels, convert pixel coordinates to centimeters, derive speed, distance to food, and zone occupancy, and score food interaction, eating, and grooming.

one_source_maze/preprocessing_synch_photometry_individual.m
three_sources_maze/preprocessing_synch_photometry_individual_sess.m
cheeseboard/preprocessing_before_or_test.m and cheeseboard/Learinig_batch.m

Preprocessing is run once per session and saves one .mat file per session, named for that session. File names carry the experimental metadata (animal ID, nutritional state, food type, and manipulation) and are parsed by the pipeline scripts, so the naming convention should be preserved.

2. Analysis and figures. The pipeline scripts load a folder of preprocessed session files, assemble them into the mice_all cell array used throughout, apply the exclusion criteria described in the manuscript, and then call the plotting and statistics functions section by section. They are written as MATLAB cell-mode scripts and are intended to be stepped through one section at a time rather than run end to end.

one_source_maze/one_source_maze_pipeline.m
three_sources_maze/three_sources_maze_pipeline.m
cheeseboard/Cheeseboard_pipeline.m

Set the data directory at the top of the pipeline script before running it. Most analysis functions take an options struct as their final argument, which sets thresholds, smoothing windows, axis limits, and figure sizes. Defaults reproduce the values used in the manuscript; the header comment of each function documents its fields. Figures are exported as vector SVG for assembly in Illustrator via export_vector_svg.

Data format

Preprocessed sessions are stored as numeric matrices, one row per tracking frame. The column layout differs slightly between modules.

GLM analysis

The models comparing predictors of AgRP activity live in one_source_maze/functions/GLM:

glm_dff_analysis.m fits models predicting z-scored dF/F from behavioral predictors for a given nutritional state and food source, with cross-validation and a circular-shift shuffle control.
glm_dff_model_comparison.m compares the candidate model set reported in Table S2.
glm_dff_per_animal.m fits the same models within individual animals.

Predictor names follow the manuscript: spatial_distance is the distance to the food source, and temporal_distance is the smaller of time since the last food encounter and time to the next one.

Data availability

Raw photometry and behavioral video are archived on Cornell servers and are available from the corresponding author on request. Preprocessed session files sufficient to run the pipeline scripts are available on the same basis.
