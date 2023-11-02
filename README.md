# value_estimation_vs_effort_mobilization
This repository contains the scripts used for the tasks and analysis of our "Value estimation versus effort mobilization: a general dissociation between ventromedial and dorsomedial prefrontal cortex" paper. Anyone who has questions is invited to contact me at nicolas.clairis at protonmail.com

To launch the task, you just have to launch one of the following three tasks within the "task" folder:
- taskLearning75.m for the learning task
- taskMentalRP.m for the Stroop task
- taskGripRP.m for the grip task.
The script will then ask you to provide a subject number, the session number (which has to be between 0 and 8 by default, script has to be modified if you want to change that) and whether you plan to use the SR ResearchLight Eyetracker or not.

To perform the analysis, you may use the following scripts depending on the analysis you want to perform:
behavior
- (under construction)
fMRI analysis
- preprocessing_NicoC_batch.m to preprocess the fMRI images obtained (some adapting of the file names and paths will be required for your own version).
- launch group_onsets_for_fMRI_MS2.m to prepare the regressors for the fMRI
- launch First_level_MS2_batch.m to obtain first level data
- launch contrasts_MS2.m to make the contrast on this first level data
- launch Second_level_batch.m to obtain Second level contrasts
- launch ConjunctionPrep_Second_level.m to perform conjunction between contrasts of interest
Pupil analysis
- preproc_eyeGazePupil_MotiScan.m for preprocessing your pupil data
- for RL task: launch MS2_pupil_RL_Val_Conf_DT.m to get correlation between Val, Conf, DT and pupil
- for grip/Stroop tasks: launch MS2_pupil_GS_X_I_RT.m to get correlation between E*, RT and pupil
