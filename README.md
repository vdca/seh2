# Deviants are detected faster at the end of verse-like sound sequences

This repository contains supplementary code and data for the following paper:

Arrazola, V. DC. (2021). Deviants are detected faster at the end of verse-like sound sequences. *Frontiers in psychology*, 12, 614872. DOI: [10.3389/fpsyg.2021.614872](https://doi.org/10.3389/fpsyg.2021.614872).

## Files

- `code/seh2_analysis.R`: script containing the functions and commands necessary to reproduce the data filtering, transformation, statistical analyses, figures and tables described in the paper.

- `data/seh2_response_data.tsv`: tab-separated dataframe containing all the reaction time (RT) responses recorded during the experiment (N=2,654), plus related variables. Each of the columns is described below.

- `data/seh2_processed_data.tsv`: tab-separated dataframe containing a subset of the observations (N=2,027); this is the dataset used in the statistical analyses. It is derived from the preceding file, but it excludes incorrect responses (e.g. missed probes, false alarms), and incorporates a number of derived variables (e.g. z-normalised RTs).

## Variables

- `subjectID`: unique identifier for each participant in the experiment.
- `itemID`: sequential identifier (1 to 76) for each trial (aka sequence) presented to a participant; it reflects the order of presentation of trials.
- `deviant`: position (1 to 8) within a sequence where a deviant tone is played (instead of a standard).
- `preceding_stds`: the number of standard tones between the current deviant and the preceding deviant; aka `deviant distance`.
- `deviant_probability`: probability that a given position contains a deviant tone (instead of a standard). Out of the 72 experimental trials presented to each participant, two thirds (=48) contain a deviant. This makes the baseline probability of any given trial of containing a deviant `p=2/3=.66`. However, this probability is lower at the beginning of a trial, and reaches `p=2/3` by the last (8th) position of the trial. Hence, the `deviant_probability` equals `p=2/3 * 1/8` on position 1, then `p=2/3 * 1/7` on position 2, and so on.
- `section`: each participant listens to 76 sequences, divided into 3 sections: training (1 to 4), block1 (5 to 40), block2 (41 to 76).
- `absRT`: the lapse of time (in ms) from the onset of a trial until the participant presses the spacebar; it's affected by a machine-specific error or lag.
- `relRT`: the lapse of time (in ms) from the onset of a deviant until the participant presses the spacebar; it's affected by a machine-specific error or lag.
- `ITI`: inter-trial interval; lapse of time (in ms) between the end of the current trial and the beginning of the following trial.
- `ISI`: inter-stimulus interval; lapse of time (in ms) between the onset of the current stimulus (=tone, stroke) and the onset of the following stimulus.
- `preISI`: interval (in ms) from the onset of the preceding stimulus to the onset of the current stimulus.
- `condition`: experimental condition (1, 2, or 3).
- `condition_label`: descriptive label for each experimental condition:
  - Cond. 1: constant ISI & ITI.
  - Cond. 2: variable ITI.
  - Cond. 3: variable ISI.
- `machine`: two experimental setups were used: Ubuntu 12.04 on a Dell XPS M1330 laptop (condition 1 and 2); Ubuntu 15.10 on Lenovo T440s laptop. The difference in hardware / software leads to systematic differences in RT lags.
- `logRT`: log-transformed version of `relRT`.
- `logRT_z`: z-normalised version of `logRT`; the RT is centered around 0, and standardised (=divided by the SD). Normalisation is performed by `machine`, because the RT lags are dependent on the setup. After z-normalising, RT values are comparable across setups.
- `relRT_z`: z-normalised version of `relRT`; similar to `logRT_z`, but based on raw RT values (in ms) instead of log-transformed RTs.
- `n_reactions`: number of times a participant presses the spacebar during a trial; the expected values are 0 times for filler trials, and 1 time for trials containing a deviant; however, participants may differ from this pattern.
- `response_type`: based on `n_reactions` and on the presence/absence of a deviant, responses are classified into four categories:
  - hit = single response to deviant;
  - overreaction = more than one response to the same deviant;
  - miss = deviant is present, but no response is recorded;
  - false alarm = a response is recorded, even if no deviant is present.
