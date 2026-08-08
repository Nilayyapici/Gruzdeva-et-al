# Gruzdeva-et-al.

This repository contains the code used to produce the analyses and figures in Gruzdeva et al., *Hunger neurons track available food locations during foraging and spatial memory recall* (2026).

The code processes fiber photometry recordings of AgRP neuron activity acquired alongside video tracking during three behavioral assays, and reproduces the quantifications, statistics, and figure panels reported in the manuscript. One complete raw session is included so that the preprocessing stage can be run without additional data.

## Requirements

- MATLAB R2019b or later
- Statistics and Machine Learning Toolbox (`fitlm`, `anova1`, `anovan`, `fitrm`, `ranova`, `multcompare`, `ranksum`, `kruskalwallis`, `ksdensity`)
- Image Processing Toolbox (`imgaussfilt`, used for spatial smoothing of occupancy and activity heatmaps)

No installation is required. Clone or download the repository, then add `common_functions`, `slanCM`, and the module you intend to run, including its `functions` subfolder, to the MATLAB path:

```matlab
addpath(genpath('common_functions'));
addpath(genpath('slanCM'));
addpath(genpath('one_source_maze'));   % or 'three_sources_maze', or 'cheeseboard'
```

## Repository organization

The repository is divided into three modules, one per behavioral assay. Each contains a preprocessing script, a top-level pipeline script, and a `functions` folder holding the analysis and plotting routines it calls. Helpers shared across modules live in `common_functions`.

| Folder | Contents |
| --- | --- |
| `one_source_maze` | Three-arm maze with a single food source; includes the GLM analyses |
| `three_sources_maze` | Three-arm maze with three sources behind mesh doors, across sessions and light/dark conditions |
| `cheeseboard` | Open cheeseboard arena, spatial learning and memory recall |
| `common_functions` | Routines used by more than one module |
| `example_data` | One complete raw session for the single-source maze |
| `slanCM` | Third-party colormap library (see Third-party code) |

`common_functions` contains `detect_food_runs`, which segments continuous tracking into individual approach runs toward and away from a food source; `findnearest`, which returns the index of the sample closest to a target time and is used throughout preprocessing to align event times to tracking frames; and `blueWhiteRed`, the diverging colormap used for the activity and difference heatmaps. `detect_food_runs` applies the same detection logic to the single-source and three-source datasets, so run-based measures are directly comparable across the two paradigms; it accepts either data layout and is configured through an optional `cfg` struct documented in its header.

## Running the code

Each module follows the same two-stage structure.

### 1. Preprocessing

The preprocessing scripts synchronize the raw acquisition streams into a single per-session matrix: photometry and analog TTL traces exported from the Doric console, centroid and DeepLabCut body-part coordinates from the top-view camera, and food-camera frame times. They compute dF/F from the 465 nm and 405 nm channels, convert pixel coordinates to centimeters, derive speed, distance to food, and zone occupancy, and incorporate the manually scored food interaction, eating, and grooming events.

- `one_source_maze/preprocessing_synch_photometry_individual.m`
- `three_sources_maze/preprocessing_synch_photometry_individual_sess.m`
- `cheeseboard/preprocessing_before_or_test.m` and `cheeseboard/Learinig_batch.m`

Preprocessing is run once per session and saves one `.mat` file per session, named for that session. File names carry the experimental metadata (animal ID, nutritional state, food type, and manipulation) and are parsed by the pipeline scripts, so the naming convention should be preserved.

To run the included example session, add example_data to the path and follow RUNNING_THE_CODE.md, which lists the file names to update before running.

```matlab
addpath('example_data');
```

The session comprises eight files:

| File | Contents |
| --- | --- |
| `example_top_cam.csv` | Frame times and centroid from the top-view camera (Bonsai, PC clock) |
| `example_DLC_top_cam.csv` | DeepLabCut labels for nose, centroid, and tail |
| `example_AI.csv` | Photometry analog input recorded in Bonsai, used for TTL synchronization |
| `example_food_cam.csv` | Frame times from the food camera |
| `photometry.mat`, `ain.mat` | Photometry and analog channels exported from the Doric console |
| `example_food_interaction.csv`, `example_eating.csv`, `example_Grooming.csv` | Manually scored event frames, as start/end pairs |

The event files hold frame indices rather than times; interaction and eating are indexed against the food camera, grooming against the top camera. The script converts both to the synchronized time base.

### 2. Analysis and figures

The pipeline scripts load a folder of preprocessed session files, assemble them into the `mice_all` cell array used throughout, apply the exclusion criteria described in the manuscript, and then call the plotting and statistics functions section by section. They are written as MATLAB cell-mode scripts and are intended to be stepped through one section at a time rather than run end to end.

- `one_source_maze/one_source_maze_pipeline.m`
- `three_sources_maze/three_sources_maze_pipeline.m`
- `cheeseboard/Cheeseboard_pipeline.m`

Set the data directory at the top of the pipeline script before running it. Most analysis functions take an `options` struct as their final argument, which sets thresholds, smoothing windows, axis limits, and figure sizes. Defaults reproduce the values used in the manuscript; the header comment of each function documents its fields. Figures are exported as vector SVG for assembly in Illustrator via `export_vector_svg`.

## Data format

Preprocessed sessions are stored as numeric matrices, one row per tracking frame. The column layout differs slightly between modules.

`one_source_maze` (11 columns):

```
1 time   2 x   3 y   4 speed   5 distance to food (shortest path)
6 path length since last visit   7 door (0 closed / 1 open)
8 food visit   9 eating (0/1)   10 grooming (0/1)   11 dF/F
```

`three_sources_maze` (13 columns):

```
1 time   2 x   3 y   4 465 nm   5 405 nm   6 dF/F   7 speed   8 zone
9 distance to food   10 distance to alternative source 2
11 distance to alternative source 1
12 door (0 closed / 1 open / 2 closed again)   13 grooming
```

`cheeseboard` (10 columns):

```
1 time   2 x   3 y   4 465 nm   5 405 nm   6 dF/F   7 speed
8 zone (0 outside, 1 food, 2 area2, 3 area3)   9 distance   10 grooming
```

## GLM analysis

The models comparing predictors of AgRP activity live in `one_source_maze/functions/GLM`:

- `glm_dff_analysis.m` fits models predicting z-scored dF/F from behavioral predictors for a given nutritional state and food source, with cross-validation and a circular-shift shuffle control.
- `glm_dff_model_comparison.m` compares the candidate model set reported in Table S2.
- `glm_dff_per_animal.m` fits the same models within individual animals.

Predictor names follow the manuscript: `spatial_distance` is the distance to the food source, and `temporal_distance` is the smaller of time since the last food encounter and time to the next one.

## Data availability

One complete raw session is included in `example_data`, sufficient to run the single-source maze preprocessing script end to end. The remaining raw photometry and behavioral video, and the full set of preprocessed session files needed to reproduce the figures, are archived on Cornell servers and are available from the corresponding author on request.

## Third-party code

`slanCM` is a colormap library by Zhaoxu Liu / slandarer, distributed through the MATLAB File Exchange and redistributed here under its own license. It is used for the `viola` and `gem` colormaps in the cheeseboard and three-source heatmaps. All other code in this repository is released under the MIT License; see `LICENSE`.

## Citation

Gruzdeva A, Shea J, Shi M, Fernandez-Ruiz A, Oliva A, Yapici N. Hunger neurons track available food locations during foraging and spatial memory recall. 2026.
