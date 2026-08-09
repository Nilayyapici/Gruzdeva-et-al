# Running the code

The scripts in this repository are session-specific templates. Each preprocessing script was written to process a single recording session, and the data folder, file names, event frame numbers, and output name for that session are written directly into the script rather than passed in as arguments. The pipeline scripts likewise carry the data folder as a literal path.

**Before running anything, set the paths and per-session values listed below.** Every script that needs editing is listed here, with the line numbers as of the current version.

Use `Figures_code.xlsx` to work backwards from a figure panel to the script and section that produce it.

---

## Before you start

Add the shared folders and the module you are working in to the MATLAB path:

```matlab
addpath(genpath('common_functions'));
addpath(genpath('slanCM'));
addpath(genpath('one_source_maze'));   % or another module
```

Every script begins with a `cd` or a `folderPath` assignment pointing at the original data location on the author's machine. Replace it with the path to your own copy before running. These are noted individually below.

---

## 1. Single-source maze

Two versions of the preprocessing script are provided.

### Example session

`one_source_maze/preprocessing_synch_photometry_individual_example.m` is set up for the recording in `one_source_maze/example_data`. All seven file names are already filled in. Only one edit is needed:

| Line | Change |
| --- | --- |
| 4 | Replace the `cd` path with the path to your copy of `one_source_maze/example_data` |

Then step through the script from the top. This is the quickest way to confirm the pipeline runs before working with your own recordings.

### Your own sessions

`one_source_maze/preprocessing_synch_photometry_individual.m` is the blank template. The `readtable` calls contain empty file names to be filled in.

| Line | Variable | Supply |
| --- | --- | --- |
| 4 | working directory | Folder holding this session's files |
| 7 | `tracking` | Top-camera frame times and centroid (Bonsai) |
| 8 | `Analog_table` | Photometry analog input (Bonsai) |
| 9 | `DLC` | DeepLabCut labels for the top camera |
| 10 | `food_cam` | Food-camera frame times |
| 211 | `interaction_food_csv` | Scored food-interaction frames |
| 214 | `eating_csv` | Scored eating frames |
| 226 | `grooming_csv` | Scored grooming frames |

`photometry.mat` and `ain.mat` (lines 15–16) are the Doric exports and keep those names for every session, so they need no editing.

Two further sets of values are session-specific:

| Line | What it is |
| --- | --- |
| 234–236 | `door_remov1`, `door1`, `door_remov2` — top-camera frame numbers at which the door was removed and replaced, read off the video. The values in the file are from the example session and are wrong for any other. |

### Scored-event file format

The three event files are comma-separated start/end frame pairs, one event per row, with no header. Food interaction and eating are indexed against the **food camera**; grooming is indexed against the **top camera**. The script converts both to the synchronized time base, so getting the two reference cameras the right way round matters.

---

## 2. Three-source maze

`three_sources_maze/preprocessing_synch_photometry_individual_sess.m`

| Line | Variable | Supply |
| --- | --- | --- |
| 5 | `tracking` | Top-camera frame times and centroid |
| 6 | `Analog_table` | Photometry analog input |
| 8 | `DLC` | DeepLabCut labels |
| 10 | `food_cam` | Food-camera frame times |
| 552 | `grooming_frames` | Grooming frames, pasted directly into the script as a numeric array rather than read from a file. Replace with the scored frames for this session. The `readtable` call on line 551 is commented out; uncomment it to supply a CSV instead, in the start/end pair format described above. 

Door open and close events should be filled in manually from the session video

---

## 3. Cheeseboard

**`cheeseboard/preprocessing_before_or_test.m`** processes one pre or test session:

| Line | Variable | Supply |
| --- | --- | --- |
| 3 | working directory | Folder holding this session's files |
| 7 | `tracking` | Top-camera frame times and centroid |
| 8 | `Analog_table` | Photometry analog input |
| 9 | `DLC` | DeepLabCut labels |
| 369 | working directory | Folder holding `Cheeseboard.mat` |
| 373–375 | index into `Cheeseboard` | Row for the animal being added. Sessions are appended to a shared cell array, so this index must be set to the correct row, or an existing session will be silently overwritten. |

**`cheeseboard/Cheeseboard_learning.m`** runs the learning analysis across animals. Set `mouseNames` and `baseDir` near line 24 to the animals to include and the top-level cheeseboard data folder. The batch-conversion steps at the top of the file are commented out and are only needed when processing new recordings; each carries its own `baseDir` to set.

Several cheeseboard functions also carry the original data path as a default argument — `analyzeAllMice_learning_cheeseboard.m`, `plotMouseTracksWithDFF.m`, `processBatchTrials.m`, `organizeFilesByMouseAndTrial.m`, and `validate_path_lengths.m`. Pass `baseDir` explicitly when calling them rather than relying on the default.

---

## 4. Object and water control

`object_water/object_water_pipeline.m` has no separate preprocessing script-the preprocessing was done using `preprocessing_synch_photometry_individual.m`; it reads session matrices produced for the single-source maze. Set `folderPath` on line 3 to the folder holding those `.mat` files.

The `sideMap` assignment near the top records which side each animal's water source was on. Extend it if you add animals.

---

## 5. Session naming convention

The pipeline scripts parse experimental metadata out of the saved file name, so the name given at the saving step determines how a session is grouped in every downstream analysis. Preserve the existing convention:

```
<animalID>_GLM_<state>_<source>.mat      e.g. F13_GLM_fasted_gel.mat
<animalID>_sess<n>.mat                    e.g. MDRE8_sess1.mat
```

The pipelines match on substrings — `fasted` / `fed` for nutritional state, `food` / `gel` for source type, and `CNO`, `opto`, `control`, `dark` for manipulation and lighting condition. A session whose file name omits one of these will be assigned to the wrong group or dropped. Save the matrix into a variable of the same name as the file, since the pipelines load each `.mat` and look for a variable matching the file name.

---

## 6. Analysis pipelines

Each pipeline needs its data directory set at the top:

- `one_source_maze/one_source_maze_pipeline.m`
- `three_sources_maze/three_sources_maze_pipeline.m` (line 3)
- `cheeseboard/Cheeseboard_pipeline.m` (line 3)
- `object_water/object_water_pipeline.m` (line 3)

They are cell-mode scripts and are meant to be stepped through section by section, not run end to end. Check the pixel-to-centimeter scale near the top of each — it is set per apparatus and per camera position and needs to be remeasured for any new setup.

---

## Checklist before each run

- [ ] `common_functions`, `slanCM`, and the module folder are on the path
- [ ] The `cd` or `folderPath` at the top of the script points at your data
- [ ] All `readtable` file names point at this session's files
- [ ] Door frame numbers match this session's video
- [ ] Grooming, eating, and interaction frames are this session's
- [ ] Pixel-to-centimeter scale matches the apparatus and camera position
- [ ] Output variable and file name follow the naming convention
