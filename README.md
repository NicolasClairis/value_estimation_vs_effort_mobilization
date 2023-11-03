# value_estimation_vs_effort_mobilization
This repository contains the scripts used for the tasks and analysis of our "Value estimation versus effort mobilization: a general dissociation between ventromedial and dorsomedial prefrontal cortex" paper. Anyone who has questions is invited to contact me at nicolas.clairis at protonmail.com

To launch the task, you just have to launch one of the following three tasks within the "task" folder:
- taskLearning75.m for the learning task
- taskMentalRP.m for the Stroop task
- taskGripRP.m for the grip task.
The script will then ask you to provide a subject number, the session number (which has to be between 0 and 8 by default, script has to be modified if you want to change that) and whether you plan to use the SR ResearchLight Eyetracker or not.

To perform the analysis, you may use the following scripts depending on the analysis you want to perform:
*behavior
- RL task:
  - MS2_RL_launch_Qmodel to perform the computational model on choices
  - MS2_RL_check_balanced_accuracy and MS2_RL_choice_model_quality to look at model quality
  - MS2_RL_choice_f_time_and_f_dQ to look at choices in function of trial number and Qleft-Qright
  - MS2_RL_avg_choices to obtain total % of "correct" choices in each condition
  - MS2_RL_Conf_DT_correl for correlation between confidence (Conf) and deliberation time (DT)
  - MS2_RL_Conf_QcQuc to check correlation between confidence and Qchosen-Qunchosen
  - MS2_RL_RT_ter: to extract graphs of reaction times (RT) in function of the variables of interest
  - MS2_RL_Val_and_Conf check how Val and Conf correlate with other variables of interest with overlap of Val and Conf on the same graphs
  - MS2_RL_Val_Conf_correl: check correlation between Val and Conf
  - MS2_RL_Val_DT_correl: check at correlation between value (Val) and deliberation time (DT)
- GS tasks:
  - MS2_GS_group_model_EV to compute the computational model on performance
  - MS2_GS_RT_f_inc to check reaction times (RT) in function of incentives (inc)
  - MS2_GS_perf_fig to check performance = f(trial) and f(incentives)


*fMRI analysis
- preprocessing_NicoC_batch.m to preprocess the fMRI images obtained (some adapting of the file names and paths will be required for your own version).
- group_onsets_for_fMRI_MS2.m to prepare the regressors for the fMRI
- First_level_MS2_batch.m to obtain first level data
- contrasts_MS2.m to make the contrast on this first level data
- Second_level_batch.m to obtain Second level contrasts
- ConjunctionPrep_Second_level.m to perform conjunction between contrasts of interest
- MS2_extract_ROI_con: script to extract ROI
- MS2_GLM1_E_RT, MS2_GS_GLM1_2_3_4, MS2_GS_GLM1_2_3_4, MS2_RL_GLM1, MS2_RL_GLM1_4_5_vmPFC_Val, MS2_RL_GLM2_3_Val_Conf_DT_orthog: script to produce the ROI figures of our manuscript once the ROI have been extracted
*Pupil analysis
- preproc_eyeGazePupil_MotiScan.m for preprocessing your pupil data
- for RL task: launch MS2_pupil_RL_Val_Conf_DT.m to get correlation between Val, Conf, DT and pupil
- for grip/Stroop tasks: launch MS2_pupil_GS_X_I_RT.m to get correlation between E*, RT and pupil
