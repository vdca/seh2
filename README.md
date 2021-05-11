# Deviants are detected faster at the end of verse-like sound sequences

This repository contains supplementary code and data for the following paper:

Arrazola, V. DC. Deviants are detected faster at the end of verse-like sound sequences.

## Files

- `code/seh2_analysis.R`: script containing the functions and commands necessary to reproduce the data filtering, transformation, statistical analyses, figures and tables described in the paper.

- `data/seh2_response_data.tsv`: tab-separated dataframe containing all the reaction time (RT) responses recorded during the experiment (N=2,654), plus related variables. Each of the columns is described below.

- `data/seh2_processed_data.tsv`: tab-separated dataframe containing a subset of the observations (N=2,027); this is the dataset used in the statistical analyses. It is derived from the preceding file, but it excludes incorrect responses (e.g. missed probes, false alarms), and incorporates a number of derived variables (e.g. z-normalised RTs).

## Variables

- `subjectID`: unique identifier for each participant in the experiment.
- `itemID`: sequential identifier (1 to 76) for each trial (aka sequence) presented to a participant; it reflects the order of presentation of trials.
- `deviant`: position (1 to 8) within a sequence where a deviant tone is played (instead of a standard).
- `section`: each participant listens to 76 sequences, divided into 3 sections: training (1 to 4), block1 (5 to 40), block2 (41 to 76).
- `absRT`: the lapse of time (in ms) from the onset of a trial until the participant presses the spacebar; it's affected by a machine-specific error or lag.
- `relRT`: the lapse of time (in ms) from the onset of a deviant until the participant presses the spacebar; it's affected by a machine-specific error or lag.
- `preceding_stds`: the number of standard tones between the current deviant and the preceding deviant; aka `deviant distance`.
- `ITI`: inter-trial interval; lapse of time (in ms) between the end of the current trial and the beginning of the following trial.
- `ISI`: inter-stimulus interval; lapse of time (in ms) between the onset of the current stimulus (=tone, stroke) and the onset of the following stimulus.
- `preISI`: interval (in ms) from the onset of the preceding stimulus to the onset of the current stimulus.
- `condition`: experimental condition (1, 2, or 3).
- `condition_label`: descriptive label for each experimental condition:
  - Cond. 1: constant ISI & ITI.
  - Cond. 2: variable ITI.
  - Cond. 3: variable ISI.
- `machine`
- `n_reactions`
- `response_type`
- `logRT`
- `logRT_z`
- `relRT_z`
- `deviant_probability`
