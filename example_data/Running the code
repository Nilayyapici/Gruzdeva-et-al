# Running the code

The preprocessing scripts in this repository are session-specific templates. Each one was written to process a single recording session, and the file names, event frame numbers, and output variable name for that session are written directly into the script rather than passed in as arguments.

**Before running any preprocessing script, edit the values listed below to match the session you are processing.** This applies to the example session included in `example_data` as well: the script as committed still names the original session files, so the file names have to be changed even to run the example.

Nothing in the analysis pipelines needs this treatment. Once sessions are preprocessed and saved, the pipeline scripts read whatever `.mat` files are in the data directory.

---

## 1. Single-source maze

`one_source_maze/preprocessing_synch_photometry_individual.m`

### File names to update

Set the current folder to the directory holding the session files, or add that folder to the MATLAB path, then replace the file names on these lines:

| Line | Variable | Replace with | For the example session |
| --- | --- | --- | --- |
| 13 | `tracking` | Top-camera frame times and centroid (Bonsai) | `example_top_cam.csv` |
| 14 | `Analog_table` | Photometry analog input (Bonsai) | `example_AI.csv` |
| 16 | `DLC` | DeepLabCut labels for the top camera | `example_DLC_top_cam.csv` |
| 18 | `food_cam` | Food-camera frame times | `example_food_cam.csv` |
| 233 | `interaction_food_csv` | Scored food-interaction frames | `example_food_interaction.csv` |
| 236 | `eating_csv` | Scored eating frames | `example_eating.csv` |
| 248 | `grooming_csv` | Scored grooming frames | `example_Grooming.csv` |

`photometry.mat` and `ain.mat` (lines 23–24) are the Doric exports and keep those names for every session, so they need no editing.

### Values to update

| Line | What it is | Note |
| --- | --- | --- |
| 256–258 | `door_remov1`, `door1`, `door_remov2` | Top-camera frame numbers at which the door was removed and replaced, read off the video. These are hardcoded for the original session and are wrong for any other. |
| 875–876 | Output variable name and save file | Uncomment and rename to follow the session-naming convention (see below). |

### Scored-event file format

The three event files are comma-separated start/end frame pairs, one event per row, with no header. Food interaction and eating are indexed against the **food camera**; grooming is indexed against the **top camera**. The script converts both to the synchronized time base, so getting the two reference cameras the right way round matters.

---

## 2. Three-source maze

`three_sources_maze/preprocessing_synch_photometry_individual_sess.m`

| Line | Variable | Replace with |
| --- | --- | --- |
| 5 | `tracking` | Top-camera frame times and centroid |
| 6 | `Analog_table` | Photometry analog input |
| 8 | `DLC` | DeepLabCut labels |
| 10 | `food_cam` | Food-camera frame times |

Additional per-session values:

- **Line 552** — grooming frames are pasted directly into the script as a numeric array rather than read from a file. Replace the array with the scored frames for the session. The `readtable` call on the line above is commented out; you can uncomment it and supply a CSV instead, in which case the format is the same start/end pair layout described above.
- **Lines 392–409** — door open and close events are derived from frame numbers set earlier in the script; check these against the session video.
- **Lines 641–648** — the saving block is commented out. Uncomment it, set the output variable name and file name, and set the destination folder.

---

## 3. Cheeseboard

Two entry points, used for different parts of the analysis.

**`cheeseboard/preprocessing_before_or_test.m`** processes one pre or test session:

| Line | Variable | Replace with |
| --- | --- | --- |
| 3 | working directory | Folder holding this session's files |
| 7 | `tracking` | Top-camera frame times and centroid |
| 8 | `Analog_table` | Photometry analog input |
| 9 | `DLC` | DeepLabCut labels |
| 369 | working directory | Folder holding `Cheeseboard.mat` |
| 373–375 | index into `Cheeseboard` | Row for the animal being added; sessions are appended to a shared cell array, so this index must be set to the correct row or an existing session will be overwritten |

**`cheeseboard/Learinig_batch.m`** runs the learning analysis across animals. Set `baseDir` on line 24 to the top-level cheeseboard data folder and `mouseNames` to the animals to include.

Several cheeseboard functions also carry the original data path as a default argument — `analyzeAllMice_learning_cheeseboard.m`, `plotMouseTracksWithDFF.m`, `processBatchTrials.m`, `organizeFilesByMouseAndTrial.m`, and `validate_path_lengths.m`. Pass `baseDir` explicitly when calling them rather than relying on the default.

---

## 4. Session naming convention

The pipeline scripts parse experimental metadata out of the saved file name, so the name given at the saving step determines how a session is grouped in every downstream analysis. Preserve the existing convention:

```
<animalID>_GLM_<state>_<source>.mat      e.g. F13_GLM_fasted_gel.mat
<animalID>_sess<n>.mat                    e.g. MDRE8_sess1.mat
```

The pipelines match on substrings — `fasted` / `fed` for nutritional state, `food` / `gel` for source type, and `CNO`, `opto`, `control`, `dark` for manipulation and lighting condition. A session whose file name omits one of these will be assigned to the wrong group or dropped. Save the matrix into a variable of the same name as the file, since the pipelines load each `.mat` and look for a variable matching the file name.

---

## 5. Analysis pipelines

The pipeline scripts need only the data directory set at the top:

- `one_source_maze/one_source_maze_pipeline.m`
- `three_sources_maze/three_sources_maze_pipeline.m` (line 3)
- `cheeseboard/Cheeseboard_pipeline.m` (line 3)

They are cell-mode scripts and are meant to be stepped through section by section, not run end to end. Check the pixel-to-centimeter scale near the top of each — it is set per apparatus and per camera position (`pix_to_cm = 0.17` for the maze, `0.16` for the cheeseboard) and needs to be remeasured if the camera has moved.

---

## Checklist before each run

- [ ] Data folder is the current directory, or on the MATLAB path
- [ ] `common_functions` and `slanCM` are on the path
- [ ] All `readtable` file names in the script point at this session's files
- [ ] Door frame numbers match this session's video
- [ ] Grooming, eating, and interaction frames are this session's
- [ ] Pixel-to-centimeter scale matches the apparatus and camera position
- [ ] Output variable and file name follow the naming convention
