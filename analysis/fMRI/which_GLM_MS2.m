function [ GLMprm, n_prm ] = which_GLM_MS2(GLM)
%[ GLMprm, n_prm ] = which_GLM_MS2(GLM)
% which_GLM_MS2 creates structure GLMprm with all the parameters of
% interest that define the GLM at the first level, contrast level and 2nd
% level of the analysis.
%
%% INPUT
% GLM: which GLM to use? (number)
%
%% OUTPUT
% GLMprm: structure containing all GLM parameters to be used in the other
% scripts (for first, second level and contrasts)
%
% (general): general condition (valid for all tasks)
% (GS) condition for both grip and stroo
% (Grip) condition concerning the grip task only
% (Stroop) condition concerning the stroop task only
% (RL) condition concerning the reinforcement learning task only
%
% n_prm: structure containing all the infos about number of parameters per
% run per task and across all
%
%% general conditions
%
%% gal.grey_mask: use grey matter mask designed during segmentation?
% (0) no: include all voxels in the analysis
%
% (1) use the individual grey mask extracted from the preprocessing in the
% first level to avoid any white matter or outside-the-brain voxels
% Check binarize_grey_mask.m and preprocessing_NicoC_batch.m for more details.
%
% (2) use a grey mask averaged across all the subjects of the study with a
% probability > X% (X depends on parameters you use in MS2_group_binarize_grey_mask.m)
% of being in the grey matter. Check MS2_group_binarize_grey_mask.m and binarize_grey_mask.m for the details.
%
% (3) use a grey mask based on SPM template for all the subjects.
% Probability of being in the grey matter at 10%.
%
% (4) use a grey mask based on SPM template for 305 subjects for all the
% subjects. Probability of being in the grey matter at 10%
%
% (5) use a mask excluding all but grey matter, based on SPM template for
% all subjects.
%
%% gal.add_drv: add temporal and/or spatial derivative
% (0) no derivative at all
% (1) time derivative added
% (2) time and dispersion derivatives added
%
%%  gal.orth_vars: orthogonalise the variables of the model or not?
% (0) avoid the orthogonalization
% (1) default option: orthogonalize the variables
%
%%  gal.zscore: zscore variables
% (0) raw values by default
% (1) zscore all variables in all tasks per run
%
%% gal.FIR:
% (1) use FIR model to model trials
%
%% gal.FIR_dur
% (X) X seconds to check after the onset
%
%% gal.FIR_nBins
% (X) the model will extract X datapoints between the onset and
% (onset+FIR_dur)
%
%% gal.onsets_only
% (1) 1 regressor/trial for each task (no modulators)
% Note: duration defined for each task in the corresponding subfield
%
%% GS conditions: will apply same parameter for both grip and stroop
% if left empty => you need/can define the parameter in the G/S
% subfield
%
%% GS.o_inc
% (1) add an onset when the incentive appears on the screen
% (2) separate gain trials from loss trials into two different regressors
% (3) separate error trials from no-error trials into two different
% regressors (pooling gain and loss conditions)
%
%% GS.dur_inc
% (0) stick
% (1) boxcar (from onset of incentive to display of effort scale)
% (2) boxcar (from onset of incentive to feedback)
%
%% GS.mod_inc.GL_cond
% (1) 0/1 for Loss/Gain condition
% (2) -1/1 for Loss/Gain condition
%
%% GS.mod_inc.R_type
% binary variable indicating whether reward presented in the form of a
% money bill or of pieces of money
% (1) 0/1 for money piece/money bill
%
%% GS.mod_inc.inc
% (1) incentive rank regressor (-6 to 6) (nominal value)
% (2) incentive value regressor (-20€ to 20€) (nominal value)
% (3) incentive rank regressor (1 to 6) (motivational value)
% (4) incentive value regressor (0.01€ to 20€) (motivational value)
%
%% GS.mod_inc.inc_bis
%(to orthogonalize to the previous one in order to get both nominal and
% motivational value in the same GLM)
% (1) incentive rank regressor (-6 to 6) (nominal value)
% (2) incentive value regressor (-20€ to 20€) (nominal value)
% (3) incentive rank regressor (1 to 6) (motivational value)
% (4) incentive value regressor (0.01€ to 20€) (motivational value)
%
%% GS.mod_inc.conf: confidence proxy
% (1) use (incentive rank - mean rank (3.5))²
% (2) use (incentive money - mean money (4.4517))²
%
%% GS.mod_inc.ROI_activity_yn
% (1) use the activity of one ROI as a regressor
%
%% GS.mod_inc.ROI_activity_GLM
% string with number of the GLM from which to extract the ROI activity ('13'/'22'/...)
%
%% GS.mod_inc.ROI_activity_period
% period when to extract ROI activity (inc/dispE/fbk)
%
%% GS.mod_inc.ROI_activity_ROI_nm
% name of the ROI to use
%
%% Grip.mod_inc.Effort
% (1) Fmax performed(i)/Fmax(i) where i is the trial number
% (2) integral of the force (normalized by Fmax(i) where i is the trial number)
% (3) Fmax performed(i)/run Fmax
% (4) integral of the force (normalized by run Fmax)
% (5) Fmax performed/estimated Fmax (by model of force)
% (6) integral of the force (normalized by estimated Fmax from the model of
%force)
%
%% Stroop.mod_inc.n_pairs_solved
% (1) number of pairs solved per trial
%
%% Stroop.mod_inc.n_incong
% (1) number of incongruent pairs solved per trial
% (2) percentage of incongruent pairs solved per trial
%
%% Stroop.mod_inc.n_errors
% (1) binary variable trial containing error(s) or not (0/1)
% (2) number of errors made per trial (only for error trials if o_inc = 3)
% (3) percentage of errors made per trial (only for error trials if o_inc = 3)
%
%% GS.mod_inc.mdl_n
% model number you want to use for extracting the cost and/or benefit
% and/or EV term and/or perf(See MS2_GS_Festimation_model_space.m for the details of
% each model)
%
%% GS.mod_inc.mdl_cost
% (1) expected cost as derived by the model of effort with actual
% performance as input
% (2) same as (1) but based on predicted performance
% (3) cost term before scaling with X/performance
%
%% GS.mod_inc.mdl_benefit
% (1) expected benefit as derived by the model of effort with actual
% performance as input
% (2) same as (1) but based on predicted performance
% (3) benefit term before scaling with performance
% (4) expected monetary payoff Val = (I*(P*G+(1-P)*L)
%
%% GS.mod_inc.mdl_EV
% (1) expected value (=benefit - cost) as derived by the model of effort with actual
% performance as input
% (2) same as (1) but based on predicted performance
%
%% Grip.mod_inc.RT_fp
% (1) raw RT first press (when starts squeezing the grip)
% (2) same as (1) but as first regressor
%
%% Stroop.mod_inc.RT_fp
% (1) raw RT first press (RT first pair)
% (2) same as (1) but as first regressor
%
%% Stroop.mod_inc.RT_mRT
% (1) mean of the raw RT of all the pairs solved after the first pair
% (=first pair RT not included for the mean because much slower)
% (2) mean of the raw RT of all the pairs solved after the first pair
% (=first pair RT not included for the mean because much slower) and also
% excluding the RT for all the errors
%
%% GS.mod_inc.totalGain_prev
% (1) total gain cumulated until previous trial in the current session
% (raw values)
%
%% GS.mod_inc.sumPerfPrev
% (1) sum of the performance (normalized between 0 and 1 for each trial)
% across all the previous trials (raw values)
% (2) sum of the performance (normalized between 0 and 1 for each trial)
% across all the previous trials (values normalized by trial number)
%
%% GS.mod_inc.lum
% (0) do not add luminance modulation
% (1) add luminance modulation as a last regressor (raw)
% (2) add luminance modulation as a first regressor (raw)
% (3) add luminance modulation as a last regressor (zscore per run)
% (4) add luminance modulation as a first regressor (zscore per run)
% (5) add luminance modulation as a last regressor (zscore across runs)
% (6) add luminance modulation as a first regressor (zscore across runs)
% (7) add luminance modulation as a last regressor (zscore across runs and per run)
% (8) add luminance modulation as a first regressor (zscore across runs and per run)
%
%% GS.mod_inc.trialN
% (1) add trial number as last regressor (raw)
% (2) add trial number as first regressor (raw)
% (3) add trial number as last regressor (zscore per run)
% (4) add trial number as first regressor (zscore per run)
% (5) add trial number middle (before RT but after performance variables) (raw)
% (6) add trial number middle (before RT but after performance variables) (zscore per run)
%
%% GS.mod_inc.inc_x_trialN
% (1) incentive*trial number (raw)
%
%% GS.mod_inc.E_pred
% (1) predicted ressource with actual performance as input (after EV)
% (2) predicted ressource with predicted performance as input (after EV)
% (3) same as (1) but placed just before benefit/cost/EV variables
% (4) same as (2) but placed just before benefit/cost/EV variables
% (5) same as (1) but placed as first regressor
% (6) same as (2) but placed as first regressor
%
%% GS.mod_inc.perf
% level of the bar at the end of the trial
% (1) performance (0-100%)
% (2) performance (0-1, normalized by 100)
% (3) performance (0-1) corrected for incentive, valence and trial number
% (4) performance (0-1) predicted by the model defined in GS.mod_inc.mdl_n
%
%% GS.o_dispE
% (1) onsets Gain+Loss pooled
% (2) onsets Gain/Loss separated
% (3) onsets error/no-error separated (for stroop mostly)
%
%% GS.dur_dispE
% (0) stick
% (1) boxcar from onset of effort scale until feedback
% (2) boxcar from onset of effort scale until RT
%
%% GS.mod_dispE.GL_cond
% (1) 0/1 for Loss/Gain condition
% (2) -1/1 for Loss/Gain condition
%
%% GS.mod_dispE.R_type
% binary variable indicating whether reward presented in the form of a
% money bill or of pieces of money
% (1) 0/1 for money piece/money bill
%
%% GS.mod_dispE.inc
% (1) incentive rank regressor (-6 to 6) (nominal value)
% (2) incentive value regressor (-20€ to 20€) (nominal value)
% (3) incentive rank regressor (1 to 6) (motivational value)
% (4) incentive value regressor (0.01€ to 20€) (motivational value)
%
%% GS.mod_dispE.inc_bis
%(to orthogonalize to the previous one in order to get both nominal and
% motivational value in the same GLM)
% (1) incentive rank regressor (-6 to 6) (nominal value)
% (2) incentive value regressor (-20€ to 20€) (nominal value)
% (3) incentive rank regressor (1 to 6) (motivational value)
% (4) incentive value regressor (0.01€ to 20€) (motivational value)
%
%% GS.mod_dispE.conf: confidence proxy
% (1) use (incentive rank - mean rank (3.5))²
% (2) use (incentive money - mean money (4.4517))²
%
%% GS.mod_dispE.ROI_activity_yn
% (1) use the activity of one ROI as a regressor
%
%% GS.mod_dispE.ROI_activity_GLM
% number of the GLM from which to extract the ROI activity (13/22/...)
%
%% GS.mod_dispE.ROI_activity_period
% period when to extract ROI activity (inc/dispE/fbk)
%
%% GS.mod_dispE.ROI_activity_ROI_nm
% name of the ROI to use
%
%% GS.mod_dispE.Effort
% hard to see how to use the same metric between stroop and grip task...
%
%% Grip.mod_dispE.Effort
% (1) Fmax performed(i)/Fmax(i) where i is the trial
% (2) integral of the force (normalized by Fmax(i))
% (3) Fmax performed(i)/run Fmax
% (4) integral of the force (normalized by run Fmax)
%
%% Stroop.mod_dispE.n_pairs_solved
% (1) number of pairs solved per trial
%
%% Stroop.mod_dispE.n_incong
% (1) number of incongruent pairs solved per trial
% (2) percentage of incongruent pairs solved per trial
%
%% Stroop.mod_dispE.n_dist
% (1) mean (number distance of the pairs solved)
% (2) mean (1/number distance of the pairs solved)
%
%% Stroop.mod_dispE.n_errors
% (1) binary variable trial containing error(s) or not (0/1)
% (2) number of errors made per trial (only for error trials if o_inc = 3)
% (3) percentage of errors made per trial (only for error trials if o_inc = 3)
%
%% GS.mod_dispE.mdl_cost
% (1) expected cost as derived by the model of effort with actual
% performance as input
% (2) same as (1) but based on predicted performance
% (3) cost term before scaling with X/performance
%
%% GS.mod_dispE.mdl_benefit
% (1) expected benefit as derived by the model of effort with actual
% performance as input
% (2) same as (1) but based on predicted performance
% (3) benefit term before scaling with performance
% (4) expected monetary payoff Val = (I*(P*G+(1-P)*L)
%
%% GS.mod_dispE.mdl_EV
% (1) expected value (=benefit - cost) as derived by the model of effort with actual
% performance as input
% (2) same as (1) but based on predicted performance
%
%% Grip.mod_dispE.RT_fp
% (1) raw RT first press (when starts squeezing the grip)
%
%% Stroop.mod_dispE.RT_fp
% (1) raw RT first press (RT first pair)
%
%% Stroop.mod_dispE.RT_mRT
% (1) mean of the raw RT of all the pairs solved after the first pair
% (=first pair RT not included for the mean because much slower)
%
%% GS.mod_dispE.lum
% think over this regressor because needs to take into account dynamic
% changes of luminance during the trial
% (0) do not add luminance modulation
% (1) luminance (incentive without effort scale) modulation (end regressor) (raw)
% (2) luminance (incentive without effort scale) modulation (1st regressor) (raw)
% (3) luminance (incentive without effort scale) modulation (end regressor) (zscore per run)
% (4) luminance (incentive without effort scale) modulation (1st regressor) (zscore per run)
% (5) luminance (incentive without effort scale) modulation (end regressor) (zscore across runs)
% (6) luminance (incentive without effort scale) modulation (1st regressor) (zscore across runs)
% (7) luminance (incentive without effort scale) modulation (end regressor) (zscore across runs and per run)
% (8) luminance (incentive without effort scale) modulation (1st regressor) (zscore across runs and per run)
%
%% GS.mod_dispE.trialN
% (1) add trial number as last regressor (raw)
% (2) add trial number as first regressor (raw)
% (3) add trial number as last regressor (zscore per run)
% (4) add trial number as first regressor (zscore per run)
%
%% GS.mod_dispE.inc_x_trialN
% (1) incentive*trial number (raw)
%
%% GS.mod_dispE.E_pred
% (1) predicted ressource with actual performance as input (after EV)
% (2) predicted ressource with predicted performance as input (after EV)
% (3) same as (1) but placed just before  benefit/cost/EV variables
% (4) same as (2) but placed just before benefit/cost/EV variables
% (5) same as (1) but placed as first regressor
% (6) same as (2) but placed as first regressor
%
%% GS.mod_dispE.perf
% level of the bar at the end of the trial
% (1) performance (0-100%)
% (2) performance (0-1, normalized by 100)
% (3) performance (0-1) corrected for incentive, valence and trial number
% (4) performance (0-1) predicted by the model defined in GS.mod_dispE.mdl_n
%
%% GS.o_perfE
% (0) do not model onset of effort separately from the onset of the display
% (1) model the onset of each performance of effort separately (=1
% onset/stroop pair for stroop and 1 onset/start of exerting effort for
% grip)
% (2) same but split gain/loss trials
% (3) same as (1) but split error/no-error trials (pooling gain+loss)
%
%% Stroop.o_perfE
% (0) do not model onset of effort separately from the onset of the display
% (1) model the onset of when the first answer is given for the first pair
% (modelling each pair is probably impossible given the timings
% (2) same but split gain/loss trials
% (3) same as (1) but split error/no-error trials (pooling gain+loss)
%
%% grip.dur_perfE
% (0) stick
% (1) boxcar from start of effort until feedback
% (2) boxcar from start of effort until end of effort
%
%% stroop.dur_perfE
% (0) stick
% (1) boxcar from RT of first pair until feedback
%
%% GS.mod_perfE.GL_cond
% (1) 0/1 for Loss/Gain condition
% (2) -1/1 for Loss/Gain condition
%
%% GS.mod_perfE.R_type
% binary variable indicating whether reward presented in the form of a
% money bill or of pieces of money
% (1) 0/1 for money piece/money bill
%
%% Grip.mod_perfE.Effort
% (1) Fmax performed/calibrated Fmax
% (2) integral of the force (normalized by calibrated Fmax)
%
%% Stroop.mod_perfE.n_pairs_solved
% (1) number of pairs solved per trial
%
%% Stroop.mod_perfE.n_incong
% (1) number of incongruent pairs per trial
% (2) percentage of incongruent pairs per trial
%
%% Stroop.mod_perfE.n_dist
% (1) mean (number distance of the pairs solved)
% (2) mean (1/number distance of the pairs solved)
%
%% Stroop.mod_perfE.n_errors
% (1) binary variable trial containing error(s) or not (0/1)
% (2) number of errors made per trial
% (3) percentage of errors made per trial
%
%% GS.mod_perfE.trialN
% (1) add trial number as last regressor (raw)
% (2) add trial number as first regressor (raw)
% (3) add trial number as last regressor (zscore per run)
% (4) add trial number as first regressor (zscore per run)
%
%% GS.o_fbk
% (1) onsets Gain+Loss
% (2) onsets Gain/Loss separated
% (3) onsets error/no-error trials separated (for stroop mostly)
%
%% GS.dur_fbk
% (0) stick
% (1) boxcar from feedback to fixation cross
%
%% GS.mod_fbk.GL_cond
% (1) 0/1 for Loss/Gain condition
% (2) -1/1 for Loss/Gain condition
%
%% GS.mod_fbk.fbk
% (1) feedback value for the current trial
%
%% GS.mod_fbk.totalGain
% (1) total gain/loss across all trials after the current trial has been
% done
%
%% GS.mod_fbk.mdl_n
% model number you want to use for extracting the cost and/or benefit
% and/or EV term and/or perf(See MS2_GS_Festimation_model_space.m for the details of
% each model)
%
%% GS.mod_fbk.E_pred
% (1) predicted ressource with actual performance as input (after EV)
% (2) predicted ressource with predicted performance as input (after EV)
% (3) same as (1) but placed just before  benefit/cost/EV variables
% (4) same as (2) but placed just before benefit/cost/EV variables
% (5) same as (1) but placed as first regressor
% (6) same as (2) but placed as first regressor
%
%% GS.mod_fbk.perf
% level of the bar at the end of the trial
% (1) performance (0-100%)
% (2) performance (0-1, normalized by 100)
% (3) performance (0-1) corrected for incentive, valence and trial number
% (4) performance (0-1) predicted by the model defined in GS.mod_inc.mdl_n
%
%% GS.mod_fbk.ROI_activity_yn
% (1) use the activity of one ROI as a regressor
%
%% GS.mod_fbk.ROI_activity_GLM
% number of the GLM from which to extract the ROI activity (13/22/...)
%
%% GS.mod_fbk.ROI_activity_period
% period when to extract ROI activity (inc/dispE/fbk)
%
%% GS.mod_fbk.ROI_activity_ROI_nm
% name of the ROI to use
%
%% GS.mod_fbk.lum
%
%
%% GS.mod_fbk.trialN
% (1) add trial number as last regressor (raw)
% (2) add trial number as first regressor (raw)
% (3) add trial number as last regressor (zscore per run)
% (4) add trial number as first regressor (zscore per run)
%
%% GS.o_cross
% (1) model cross onsets separately for grip and stroop task
% (2) model cross onsets pooling grip and stroop task
%
%% GS.dur_cross
% (0) stick
% (1) boxcar from cross onset to incentive display
%
%% GS.o_missed_trials_dispE
% (1) missed trials onset
%
%% GS.dur_missed_trials_dispE
% (0) stick at display effort scale onset
% (1) boxcar from onset of display effort scale to feedback
%
%
%
%% RL conditions
%% RL.o_stim
% (1) add an onset when the pair of stimuli appears on screen
% (2) separate onsets depending on pair type: gain/neutral/loss pair
% (3) pool gain and loss pair together but put neutral pair apart
% (4) separate onsets depending on trials (first/second half)
% (5) separate onsets depending on pair type (gain/neutral/loss) and trials
% (first/second half)
% (6) pool gain and loss pair together but put neutral pair apart + split
% first/second half of trials
%
%% RL.dur_stim
% (0) stick function
% (1) boxcar from stimulus display until choice in red
% (2) boxcar from stimulus display until RT
% (3) boxcar from stimulus display until feedback
%
%% RL.mod_stim.trialN
% (1) trial number modulation as first regressor (raw)
% (2) zscore trial number per run as first regressor
% (3) zscore trial number across runs as first regressor
% (4) trial number (raw) just before ROI and RT regressors
% (5) zscore trial number per run just before ROI and RT regressors
% (6) zscore trial number across runs just before ROI and RT regressors
% (7) trial number (raw) as last regressor
% (8) zscore trial number per run  as last regressor
% (9) zscore trial number across runs  as last regressor
%
%% RL.mod_stim.mdl_type
% 'Fabien': based on MS2_RL_Fabien_version.m output
% 'Nico': based on MS2_RL_model_bis.m output
% 'Nico_avg_prm': based on MS2_RL_model_bis.m average parameters for each
% subject
%
%% RL.mod_stim.mdl_n: model number to use
%
%% RL.mod_stim.SV
% (1) SV = pA*QA+pB*QB (raw values)
% (2) SV = QA+QB (raw values)
% (3) SV = Qchosen - Qunchosen
% (4) SV = Qchosen
% (5) SV = Qchosen/(QA+QB)
% (6-10): same as (1-5) in the same order, but as a last regressor (after
% Conf and DT)
%
%% RL.mod_stim.dQ (if RL.o_stim = 2 then only use this for gain and loss pair,
% but not for neutral pair since dQ = 0 always in that case)
% model 3:
% (1)  QbestCue/pair - QnonbestCue/pair (raw values based on VBA model 3)
% (2)  same as (1) but values zscored/run
% (3)  same as (1) but values zscored across runs
% (4)  same as (1) but values zscored/run/pair
% (5)  same as (1) but values zscored across runs per pair
% (6)  [QbestCue-QnonBestCue]^2 (raw values based on VBA model 3)
% (7) |QbestCue - QnonBestCue|/sqrt(QbestCue^2 + QnonBestcue^2)
% (8) same as (7) but zscored/run
%
%% RL.mod_stim.pBest
% (1) ( p(choice = best option) - 0.5)^2
% (2) ( p(choice = left option) - 0.5)^2
% (3) p(choice = chosen option)
%
%% RL.mod_stim.ROI_activity_yn
% (1) use the activity of one ROI as a regressor
%
%% RL.mod_stim.ROI_activity_GLM
% number of the GLM from which to extract the ROI activity (13/22/...)
%
%% RL.mod_stim.ROI_activity_period
% period when to extract ROI activity (stim/fbk)
%
%% RL.mod_stim.ROI_activity_ROI_nm
% name of the ROI to use
%
%% RL.mod_stim.RT
% (1) raw RT
% (2) RT zscored per RL run
% (3) RT zscored across all RL runs (but RL runs only)
% (4) RT zscored across all RL runs + per run
% (5) RT corrected for several factors
% (6) RT corrected for several factors + zscored per run
% (7) RT corrected for several factors + zscored across RL runs
% (8) RT corrected for several factors + zscored across RL runs and per run
% (9-16) same as (1-8) (same order) but as first regressor and not last
%
%% RL.o_answer
% (1) add an onset for when the choice is being done
% (2) separate choice onsets depending on pair type: gain/neutral/loss pair
%
%% RL.dur_answer
% (0) stick
% (1) boxcar from RT until choice in red
% (2) boxcar from RT until feedback
%
%% RL.o_chosen
% (1) add an onset when the chosen stimulus appears in red on screen
% (2) separate onsets depending on pair type: gain/neutral/loss pair
% (3) pool gain and loss pair together but put neutral pair apart
% (4) separate onsets depending on trials (first/second half)
% (5) separate onsets depending on pair type (gain/neutral/loss) and trials
% (first/second half)
%
%% RL.dur_chosen
% (0) stick function
% (1) boxcar from choice in red until feedback
%
%% RL.mod_chosen.trialN
% (1) trial number modulation as first regressor (raw)
% (2) zscore trial number per run as first regressor
% (3) zscore trial number across runs as first regressor
% (4) trial number (raw) just before ROI and RT regressors
% (5) zscore trial number per run just before ROI and RT regressors
% (6) zscore trial number across runs just before ROI and RT regressors
% (7) trial number (raw) as last regressor
% (8) zscore trial number per run  as last regressor
% (9) zscore trial number across runs  as last regressor
%
%% RL.mod_chosen.mdl_type
% 'Fabien': based on MS2_RL_Fabien_version.m output
% 'Nico': based on MS2_RL_model_bis.m output
% 'Nico_avg_prm': based on MS2_RL_model_bis.m average parameters for each
% subject
%
%% RL.mod_chosen.mdl_n: model number to use
%
%% RL.mod_chosen.SV
% (1) SV = pA*QA+pB*QB (raw values)
% (2) SV = QA+QB (raw values)
% (3) SV = Qchosen - Qunchosen (raw values)
% (4) SV = Qchosen (raw values)
% (5) SV = Qchosen/(QA+QB) (raw values)
%
%% RL.mod_chosen.dQ (if RL.o_chosen = 2 then only use this for gain and loss pair,
% but not for neutral pair since dQ = 0 always in that case)
% (1)  QbestCue/pair - QnonbestCue/pair (raw values based on VBA model 3)
% (2)  same as (1) but values zscored/run
% (3)  same as (1) but values zscored across runs
% (4)  same as (1) but values zscored/run/pair
% (5)  same as (1) but values zscored across runs per pair
% (6)  [QbestCue-QnonBestCue]^2 (raw values based on VBA model 3)
% (7) |QbestCue - QnonBestCue|/sqrt(QbestCue^2 + QnonBestcue^2)
% (8) same as (7) but zscored/run
%
%% RL.mod_chosen.pBest
% (1) ( p(choice = best option) - 0.5)^2
% (2) ( p(choice = left option) - 0.5)^2
% (3) p(choice = chosen option)
%
%% RL.mod_chosen.ROI_activity_yn
% (1) use the activity of one ROI as a regressor
%
%% RL.mod_chosen.ROI_activity_GLM
% number of the GLM from which to extract the ROI activity (13/22/...)
%
%% RL.mod_chosen.ROI_activity_period
% period when to extract ROI activity (stim/fbk)
%
%% RL.mod_chosen.ROI_activity_ROI_nm
% name of the ROI to use
%
%% RL.mod_chosen.RT
% (1) raw RT
% (2) RT zscored per RL run
% (3) RT zscored across all RL runs (but RL runs only)
% (4) RT zscored across all RL runs + per run
% (5) RT corrected for several factors
% (6) RT corrected for several factors + zscored per run
% (7) RT corrected for several factors + zscored across RL runs
% (8) RT corrected for several factors + zscored across RL runs and per run
%
%% RL.o_fbk
% (0) do not model onset of feedback
% (1) add an onset for all feedbacks (grouped)
% (2) separate feedback onsets depending on pair type gain/neutral/loss
% pair
% (3) separate feedback onsets depending on feedback type gain/neutral/loss
% (4) separate feedback onsets depending on pair type gain/neutral/loss and
% also for each pair, depending on feedback type gain/neutral/loss
% (5) pool gain and loss but leave neutral pair apart
% (6) pool gain and loss but leave neutral pair apart
% and separate onsets for gain and loss only depending on trials (first/second half)
%
%% RL.dur_fbk
% (0) stick
% (1) boxcar from feedback onset to cross display
%
%% RL.mod_fbk.trialN
% (1) trial number modulation
% (2) zscore trial number per run
% (3) zscore trial number across runs
%
%% RL.mod_fbk.fbk
% (1) feedback (current trial) as modulator
%
%% RL.mod_fbk.mdl_type
% 'Fabien': based on MS2_RL_Fabien_version.m output
% 'Nico': based on MS2_RL_model_bis.m output
% 'Nico_avg_prm': based on MS2_RL_model_bis.m average parameters for each
% subject
%
%% RL.mod_fbk.mdl_n: model number to use
%
%% RL.mod_fbk.PE
% (1) Prediction Error modulator (based on model defined by
% RL.mod_fbk.mdl_type and RL.mod_fbk.mdl_n)
% (2) |Prediction Error| modulator (based on model defined by
% RL.mod_fbk.mdl_type and RL.mod_fbk.mdl_n)
%
%% RL.mod_fbk.PE_bis
% (1) Prediction Error modulator (based on model defined by
% RL.mod_fbk.mdl_type and RL.mod_fbk.mdl_n)
% (2) |Prediction Error| modulator (based on model defined by
% RL.mod_fbk.mdl_type and RL.mod_fbk.mdl_n)
%
%% RL.mod_fbk.totalGain
% (1) total gain/loss across all trials after the current trial has been
% done
%
%% RL.mod_fbk.ROI_activity_yn
% (1) use the activity of one ROI as a regressor
%
%% RL.mod_fbk.ROI_activity_GLM
% number of the GLM from which to extract the ROI activity (13/22/...)
%
%% RL.mod_fbk.ROI_activity_period
% period when to extract ROI activity (stim/fbk)
%
%% RL.mod_fbk.ROI_activity_ROI_nm
% name of the ROI to use
%
%% RL.o_cross
% (1) add fixation cross onsets
%
%% RL.dur_cross
% (0) stick
% (1) boxcar from cross onset until stim display
%
%% RL.o_missed_trials_stim
% (1) add missed trials onsets
%
%% RL.dur_missed_trials_stim
% (0) stick
% (1) duration from onset to feedback
%
%
% See also binarize_grey_mask.m and preprocessing_NicoC_batch.m

%% list of all possible onsets and regressors
[list_all_GLMprm, n_potReg] = MS2_list_potential_onsets_and_modulators;

%% set all parameters to zero by default

% general
for iGal = 1:n_potReg.gal
    curr_gal_prm = list_all_GLMprm.gal{iGal};
    gal.(curr_gal_prm) = 0;
end
gal.orth_vars = 1; % default = orthogonalize variables

% RL
% onsets & durations
for iORL = 1:n_potReg.RL.onsets
    curr_RL_onset_prm   = list_all_GLMprm.RL.onsets{iORL};
    curr_RL_dur_prm     = list_all_GLMprm.RL.durations{iORL};
    RL.(curr_RL_onset_prm) = 0;
    RL.(curr_RL_dur_prm) = 0;
end
% loop through modulators
for i_RL_mod = 1:n_potReg.RL.mods
    curr_mod = list_all_GLMprm.RL.modulators{i_RL_mod};
    for iRL_regs = 1:n_potReg.RL.(curr_mod) % loop through regressors
        curr_RL_prm = list_all_GLMprm.RL.(curr_mod){iRL_regs};
        RL.(curr_mod).(curr_RL_prm) = 0;
    end
end

% special case for model type which is a string and not a numerical value
RL.mod_stim.mdl_type    = '';
RL.mod_chosen.mdl_type  = '';
RL.mod_fbk.mdl_type     = '';

% grip
% onsets
for iOGrip = 1:n_potReg.grip.onsets
    curr_grip_onset_prm = list_all_GLMprm.grip.onsets{iOGrip};
    curr_grip_dur_prm   = list_all_GLMprm.grip.durations{iOGrip};
    grip.(curr_grip_onset_prm) = 0;
    grip.(curr_grip_dur_prm) = 0;
end
% loop through modulators
for i_G_mod = 1:n_potReg.grip.mods
    curr_Gmod = list_all_GLMprm.grip.modulators{i_G_mod};
    for iG_regs = 1:n_potReg.grip.(curr_Gmod)
        curr_grip_prm = list_all_GLMprm.grip.(curr_Gmod){iG_regs};
        grip.(curr_Gmod).(curr_grip_prm) = 0;
    end
end

% stroop
% onsets
for iOStroop = 1:n_potReg.stroop.onsets
    curr_stroop_onset_prm   = list_all_GLMprm.stroop.onsets{iOStroop};
    curr_stroop_dur_prm     = list_all_GLMprm.stroop.durations{iOStroop};
    stroop.(curr_stroop_onset_prm) = 0;
    stroop.(curr_stroop_dur_prm) = 0;
end
% modulators
for i_S_mod = 1:n_potReg.stroop.mods
    curr_S_mod = list_all_GLMprm.stroop.modulators{i_S_mod};
    for iS_regs = 1:n_potReg.stroop.(curr_S_mod)
        curr_stroop_prm = list_all_GLMprm.stroop.(curr_S_mod){iS_regs};
        stroop.(curr_S_mod).(curr_stroop_prm) = 0;
    end
end

%% define GLM properties (when needs to be different than zero)
switch GLM
    case 1
        % RL: Val/Conf/DT
        % GS: check X/RT no orthogonalization
        gal.add_drv = 0;
        gal.grey_mask = 0;
        gal.orth_vars = 0;
        gal.zscore = 1;
        % grip
        grip.o_inc = 1;
        grip.dur_inc = 0;
        grip.mod_inc.mdl_n = 1;
        grip.mod_inc.E_pred = 6;
        grip.mod_inc.RT_fp = 1;
        grip.o_dispE = 1;
        grip.dur_dispE = 1;
        grip.o_fbk = 1;
        grip.dur_fbk = 0;
        % stroop
        stroop.o_inc = 1;
        stroop.dur_inc = 0;
        stroop.mod_inc.mdl_n = 1;
        stroop.mod_inc.E_pred = 6;
        stroop.mod_inc.RT_fp = 1;
        stroop.o_dispE = 1;
        stroop.dur_dispE = 1;
        stroop.o_fbk = 1;
        stroop.dur_fbk = 0;
        % RL
        RL.o_stim = 3;
        RL.dur_stim = 1;
        RL.mod_stim.mdl_type = 'Nico';
        RL.mod_stim.mdl_n = 6;
        RL.mod_stim.SV = 1;
        RL.mod_stim.pBest = 2;
        RL.mod_stim.RT = 1;
        RL.o_chosen = 1;
        RL.dur_chosen = 1;
        RL.o_fbk = 5;
        RL.dur_fbk = 0;
        RL.mod_fbk.mdl_type = RL.mod_stim.mdl_type;
        RL.mod_fbk.mdl_n = RL.mod_stim.mdl_n; % same as for stim
        RL.mod_fbk.PE = 1;
        %% series of control GLMs for GLM 1
    case 2
        % regressors orthogonalized in the following order:
        % RL: Val/Conf/DT
        % GS: X/RT orthogonalized
        gal.add_drv = 0;
        gal.grey_mask = 0;
        gal.orth_vars = 1;
        gal.zscore = 1;
        % grip
        grip.o_inc = 1;
        grip.dur_inc = 0;
        grip.mod_inc.mdl_n = 1;
        grip.mod_inc.E_pred = 6;
        grip.mod_inc.RT_fp = 1;
        grip.o_dispE = 1;
        grip.dur_dispE = 1;
        grip.o_fbk = 1;
        grip.dur_fbk = 0;
        % stroop
        stroop.o_inc = 1;
        stroop.dur_inc = 0;
        stroop.mod_inc.mdl_n = 1;
        stroop.mod_inc.E_pred = 6;
        stroop.mod_inc.RT_fp = 1;
        stroop.o_dispE = 1;
        stroop.dur_dispE = 1;
        stroop.o_fbk = 1;
        stroop.dur_fbk = 0;
        % RL
        RL.o_stim = 3;
        RL.dur_stim = 1;
        RL.mod_stim.mdl_type = 'Nico';
        RL.mod_stim.mdl_n = 6;
        RL.mod_stim.SV = 1;
        RL.mod_stim.pBest = 2;
        RL.mod_stim.RT = 1;
        RL.o_chosen = 1;
        RL.dur_chosen = 1;
        RL.o_fbk = 5;
        RL.dur_fbk = 0;
        RL.mod_fbk.mdl_type = RL.mod_stim.mdl_type;
        RL.mod_fbk.mdl_n = RL.mod_stim.mdl_n; % same as for stim
        RL.mod_fbk.PE = 1;
    case 3
        % reverted serial orthogonalization (like GLM 2 but reverted order)
        % RL: DT/Conf/Val
        % GS: RT/X
        gal.add_drv = 0;
        gal.grey_mask = 0;
        gal.orth_vars = 1;
        gal.zscore = 1;
        % grip
        grip.o_inc = 1;
        grip.dur_inc = 0;
        grip.mod_inc.mdl_n = 1;
        grip.mod_inc.E_pred = 2;
        grip.mod_inc.RT_fp = 2;
        grip.o_dispE = 1;
        grip.dur_dispE = 1;
        grip.o_fbk = 1;
        grip.dur_fbk = 0;
        % stroop
        stroop.o_inc = 1;
        stroop.dur_inc = 0;
        stroop.mod_inc.mdl_n = 1;
        stroop.mod_inc.E_pred = 2;
        stroop.mod_inc.RT_fp = 2;
        stroop.o_dispE = 1;
        stroop.dur_dispE = 1;
        stroop.o_fbk = 1;
        stroop.dur_fbk = 0;
        % RL
        RL.o_stim = 3;
        RL.dur_stim = 1;
        RL.mod_stim.mdl_type = 'Nico';
        RL.mod_stim.mdl_n = 6;
        RL.mod_stim.SV = 6;
        RL.mod_stim.pBest = 2;
        RL.mod_stim.RT = 9;
        RL.o_chosen = 1;
        RL.dur_chosen = 1;
        RL.o_fbk = 5;
        RL.dur_fbk = 0;
        RL.mod_fbk.mdl_type = RL.mod_stim.mdl_type;
        RL.mod_fbk.mdl_n = RL.mod_stim.mdl_n; % same as for stim
        RL.mod_fbk.PE = 1;
    
    case 4
        % RL: use Qa+Qb for Val
        % GS: use X/RT during performance instead of incentive
        gal.add_drv = 0;
        gal.grey_mask = 0;
        gal.orth_vars = 0;
        gal.zscore = 1;
        % grip
        grip.o_inc = 1;
        grip.dur_inc = 0;
        grip.o_dispE = 1;
        grip.dur_dispE = 1;
        grip.mod_dispE.mdl_n = 1;
        grip.mod_dispE.E_pred = 6;
        grip.mod_dispE.RT_fp = 1;
        grip.o_fbk = 1;
        grip.dur_fbk = 0;
        % stroop
        stroop.o_inc = 1;
        stroop.dur_inc = 0;
        stroop.o_dispE = 1;
        stroop.dur_dispE = 1;
        stroop.mod_dispE.mdl_n = 1;
        stroop.mod_dispE.E_pred = 6;
        stroop.mod_dispE.RT_fp = 1;
        stroop.o_fbk = 1;
        stroop.dur_fbk = 0;
        % RL
        RL.o_stim = 3;
        RL.dur_stim = 1;
        RL.mod_stim.mdl_type = 'Nico';
        RL.mod_stim.mdl_n = 6;
        RL.mod_stim.SV = 2;
        RL.mod_stim.pBest = 2;
        RL.mod_stim.RT = 1;
        RL.o_chosen = 1;
        RL.dur_chosen = 1;
        RL.o_fbk = 5;
        RL.dur_fbk = 0;
        RL.mod_fbk.mdl_type = RL.mod_stim.mdl_type;
        RL.mod_fbk.mdl_n = RL.mod_stim.mdl_n; % same as for stim
        RL.mod_fbk.PE = 1;
    case 5
        % RL: use Qch-Quc for Val
        gal.add_drv = 0;
        gal.grey_mask = 0;
        gal.orth_vars = 0;
        gal.zscore = 1;
        % grip
        grip.o_inc = 1;
        grip.dur_inc = 0;
        grip.mod_inc.mdl_n = 1;
        grip.mod_inc.E_pred = 6;
        grip.mod_inc.RT_fp = 1;
        grip.o_dispE = 1;
        grip.dur_dispE = 1;
        grip.o_fbk = 1;
        grip.dur_fbk = 0;
        % stroop
        stroop.o_inc = 1;
        stroop.dur_inc = 0;
        stroop.mod_inc.mdl_n = 1;
        stroop.mod_inc.E_pred = 6;
        stroop.mod_inc.RT_fp = 1;
        stroop.o_dispE = 1;
        stroop.dur_dispE = 1;
        stroop.o_fbk = 1;
        stroop.dur_fbk = 0;
        % RL
        RL.o_stim = 3;
        RL.dur_stim = 1;
        RL.mod_stim.mdl_type = 'Nico';
        RL.mod_stim.mdl_n = 6;
        RL.mod_stim.SV = 3;
        RL.mod_stim.pBest = 2;
        RL.mod_stim.RT = 1;
        RL.o_chosen = 1;
        RL.dur_chosen = 1;
        RL.o_fbk = 5;
        RL.dur_fbk = 0;
        RL.mod_fbk.mdl_type = RL.mod_stim.mdl_type;
        RL.mod_fbk.mdl_n = RL.mod_stim.mdl_n; % same as for stim
        RL.mod_fbk.PE = 1;
        
        %% control GLM for GLM 1: competition replace E* by |incentives| to compare betas later
    case 6 % check |Inc|/RT no orthogonalization
        gal.add_drv = 0;
        gal.grey_mask = 0;
        gal.orth_vars = 0;
        gal.zscore = 1;
        % grip
        grip.o_inc = 1;
        grip.dur_inc = 0;
        grip.mod_inc.inc = 4;
        grip.mod_inc.RT_fp = 1;
        grip.o_dispE = 1;
        grip.dur_dispE = 1;
        grip.o_fbk = 1;
        grip.dur_fbk = 0;
        % stroop
        stroop.o_inc = 1;
        stroop.dur_inc = 0;
        stroop.mod_inc.inc = 4;
        stroop.mod_inc.RT_fp = 1;
        stroop.o_dispE = 1;
        stroop.dur_dispE = 1;
        stroop.o_fbk = 1;
        stroop.dur_fbk = 0;
        % RL
        RL.o_stim = 3;
        RL.dur_stim = 1;
        RL.mod_stim.mdl_type = 'Nico';
        RL.mod_stim.mdl_n = 6;
        RL.mod_stim.SV = 1;
        RL.mod_stim.pBest = 2;
        RL.mod_stim.RT = 1;
        RL.o_chosen = 1;
        RL.dur_chosen = 1;
        RL.o_fbk = 5;
        RL.dur_fbk = 0;
        RL.mod_fbk.mdl_type = RL.mod_stim.mdl_type;
        RL.mod_fbk.mdl_n = RL.mod_stim.mdl_n; % same as for stim
        RL.mod_fbk.PE = 1;
end

%% store in output
GLMprm.gal = gal;
GLMprm.RL = RL;

% if GS used, copy parameter inside grip and stroop fields as well
% be careful: GS to be used only if no difference at all between two
% conditions
if exist('GS','var') && ~isempty(GS)
    GLMprm.grip = GS;
    GLMprm.stroop = GS;
else
    GLMprm.grip = grip;
    GLMprm.stroop = stroop;
end

%% check incompatible parameters

if ismember(grip.mod_inc.Effort,5:6)
    error('case not ready yet, please update scripts and then do it');
end

% predicted performance and model
if grip.mod_inc.perf == 4 && grip.mod_inc.mdl_n == 0
    error(['You need to define the model number in grip.mod_inc.mdl_n ',...
        'to be able to use the predicted performance of the model']);
end
if stroop.mod_inc.perf == 4 && stroop.mod_inc.mdl_n == 0
    error(['You need to define the model number in stroop.mod_inc.mdl_n ',...
        'to be able to use the predicted performance of the model']);
end
if grip.mod_dispE.perf == 4 && grip.mod_dispE.mdl_n == 0
    error(['You need to define the model number in grip.mod_dispE.mdl_n ',...
        'to be able to use the predicted performance of the model']);
end
if grip.mod_dispE.perf == 4 && grip.mod_dispE.mdl_n == 0
    error(['You need to define the model number in stroop.mod_dispE.mdl_n ',...
        'to be able to use the predicted performance of the model']);
end
if grip.mod_fbk.perf == 4 && grip.mod_fbk.mdl_n == 0
    error(['You need to define the model number in grip.mod_fbk.mdl_n ',...
        'to be able to use the predicted performance of the model']);
end
if stroop.mod_fbk.perf == 4 && stroop.mod_fbk.mdl_n == 0
    error(['You need to define the model number in stroop.mod_fbk.mdl_n ',...
        'to be able to use the predicted performance of the model']);
end

% split stroop trials by presence or not of errors and add on top of that
% binary modulation by error presence should not be done since the
% modulator and the trial categorization will completely overlap
if stroop.o_inc == 3 && stroop.mod_inc.n_errors == 1
    error(['stroop.o_inc == 3 && stroop.mod_inc.n_errors == 1',...
        ' should not be.']);
end
if stroop.o_dispE == 3 && stroop.mod_dispE.n_errors == 1
    error(['stroop.o_perfE == 3 && stroop.mod_perfE.n_errors == 1',...
        ' should not be.']);
end
if stroop.o_perfE == 3 && stroop.mod_perfE.n_errors == 1
    error(['stroop.o_perfE == 3 && stroop.mod_perfE.n_errors == 1',...
        ' should not be.']);
end

% if GL_cond > 0, conditions should not be split
if (grip.mod_inc.GL_cond ~= 0 && grip.o_inc == 2)
    error('grip.mod_inc.GL_cond>0 and grip.o_inc = 2 should not be allowed');
end
if (stroop.mod_inc.GL_cond ~= 0 && stroop.o_inc == 2)
    error('stroop.mod_inc.GL_cond>0 and stroop.o_inc = 2 should not be allowed');
end
if (grip.mod_dispE.GL_cond ~= 0 && grip.o_dispE == 2)
    error('grip.mod_dispE.GL_cond>0 and grip.o_dispE = 2 should not be allowed');
end
if (stroop.mod_dispE.GL_cond ~= 0 && stroop.o_dispE == 2)
    error('stroop.mod_dispE.GL_cond>0 and stroop.o_dispE = 2 should not be allowed');
end
if (grip.mod_perfE.GL_cond ~= 0 && grip.o_perfE == 2)
    error('grip.mod_perfE.GL_cond>0 and grip.o_perfE = 2 should not be allowed');
end
if (stroop.mod_perfE.GL_cond ~= 0 && stroop.o_perfE == 2)
    error('stroop.mod_perfE.GL_cond>0 and stroop.o_perfE = 2 should not be allowed');
end
if (grip.mod_fbk.GL_cond ~= 0 && grip.o_fbk == 2)
    error('grip.mod_fbk.GL_cond>0 and grip.o_fbk = 2 should not be allowed');
end
if (stroop.mod_fbk.GL_cond ~= 0 && stroop.o_fbk == 2)
    error('stroop.mod_fbk.GL_cond>0 and stroop.o_fbk = 2 should not be allowed');
end

% incentive
if (grip.mod_inc.inc == grip.mod_inc.inc_bis && grip.mod_inc.inc ~= 0) ||...
        (grip.mod_dispE.inc == grip.mod_dispE.inc_bis && grip.mod_dispE.inc ~= 0) ||...
        (grip.mod_perfE.inc == grip.mod_perfE.inc_bis && grip.mod_perfE.inc ~= 0) ||...
        (stroop.mod_inc.inc == stroop.mod_inc.inc_bis && stroop.mod_inc.inc ~= 0) ||...
        (stroop.mod_dispE.inc == stroop.mod_dispE.inc_bis && stroop.mod_dispE.inc ~= 0) ||...
        (stroop.mod_perfE.inc == stroop.mod_perfE.inc_bis && stroop.mod_perfE.inc ~= 0)
    error('inc = inc_bis should not happen. Please fix that.');
end

% if onsets_only, no modulator should be allowed
if gal.onsets_only ~= 0
    
    RL_mod_stim_fields      = fieldnames(RL.mod_stim);
    for iField = 1:length(RL_mod_stim_fields)
        mod_prm = RL.mod_stim.(RL_mod_stim_fields{iField});
        if mod_prm ~= 0
            error(['onsets_only GLM but RL.mod_stim.',RL_mod_stim_fields{iField},' = ',num2str(mod_prm)]);
        end
    end
    
    RL_mod_fbk_fields       = fieldnames(RL.mod_fbk);
    for iField = 1:length(RL_mod_fbk_fields)
        mod_prm = RL.mod_fbk.(RL_mod_fbk_fields{iField});
        if mod_prm ~= 0
            error(['onsets_only GLM but RL.mod_fbk.',RL_mod_fbk_fields{iField},' = ',num2str(mod_prm)]);
        end
    end
    
    G_mod_inc_fields        = fieldnames(grip.mod_inc);
    for iField = 1:length(G_mod_inc_fields)
        mod_prm = grip.mod_inc.(G_mod_inc_fields{iField});
        if mod_prm ~= 0
            error(['onsets_only GLM but grip.mod_inc.',G_mod_inc_fields{iField},' = ',num2str(mod_prm)]);
        end
    end
    
    G_mod_dispE_fields      = fieldnames(grip.mod_dispE);
    for iField = 1:length(G_mod_dispE_fields)
        mod_prm = grip.mod_dispE.(G_mod_dispE_fields{iField});
        if mod_prm ~= 0
            error(['onsets_only GLM but grip.mod_dispE.',G_mod_dispE_fields{iField},' = ',num2str(mod_prm)]);
        end
    end
    
    G_mod_perfE_fields      = fieldnames(grip.mod_perfE);
    for iField = 1:length(G_mod_perfE_fields)
        mod_prm = grip.mod_perfE.(G_mod_perfE_fields{iField});
        if mod_prm ~= 0
            error(['onsets_only GLM but grip.mod_perfE.',G_mod_perfE_fields{iField},' = ',num2str(mod_prm)]);
        end
    end
    
    G_mod_fbk_fields        = fieldnames(grip.mod_fbk);
    for iField = 1:length(G_mod_fbk_fields)
        mod_prm = grip.mod_fbk.(G_mod_fbk_fields{iField});
        if mod_prm ~= 0
            error(['onsets_only GLM but grip.mod_fbk.',G_mod_fbk_fields{iField},' = ',num2str(mod_prm)]);
        end
    end
    
    S_mod_inc_fields        = fieldnames(stroop.mod_inc);
    for iField = 1:length(S_mod_inc_fields)
        mod_prm = stroop.mod_inc.(S_mod_inc_fields{iField});
        if mod_prm ~= 0
            error(['onsets_only GLM but stroop.mod_inc.',S_mod_inc_fields{iField},' = ',num2str(mod_prm)]);
        end
    end
    
    S_mod_dispE_fields      = fieldnames(stroop.mod_dispE);
    for iField = 1:length(S_mod_dispE_fields)
        mod_prm = stroop.mod_dispE.(S_mod_dispE_fields{iField});
        if mod_prm ~= 0
            error(['onsets_only GLM but stroop.mod_dispE.',S_mod_dispE_fields{iField},' = ',num2str(mod_prm)]);
        end
    end
    
    S_mod_perfE_fields      = fieldnames(stroop.mod_perfE);
    for iField = 1:length(S_mod_perfE_fields)
        mod_prm = stroop.mod_perfE.(S_mod_perfE_fields{iField});
        if mod_prm ~= 0
            error(['onsets_only GLM but stroop.mod_perfE.',S_mod_perfE_fields{iField},' = ',num2str(mod_prm)]);
        end
    end
    
    S_mod_fbk_fields        = fieldnames(stroop.mod_fbk);
    for iField = 1:length(S_mod_fbk_fields)
        mod_prm = stroop.mod_fbk.(S_mod_fbk_fields{iField});
        if mod_prm ~= 0
            error(['onsets_only GLM but stroop.mod_fbk.',S_mod_fbk_fields{iField},' = ',num2str(mod_prm)]);
        end
    end
end

% FIR
if gal.FIR == 0 &&...
        (gal.FIR_dur ~= 0  || gal.FIR_nBins ~= 0)
    error('gal.FIR = 0, but one of the FIR parameters in one of the tasks is not equal to zero. Please fix it.');
end
if gal.FIR ~= 0 &&...
        (gal.FIR_dur == 0  || gal.FIR_nBins == 0)
    error('gal.FIR different from zero but one of the FIR parameters is equal to zero in at least one task. Please fix it');
end

%% output GLM depending on the selected factors + number of parameters per task
which_GLM_MS2_infos_gal( GLMprm);
fprintf('\n \n');

n_prm_RL_per_run = which_GLM_MS2_infos_RL( GLMprm );
fprintf('\n \n');

n_prm_grip_per_run = which_GLM_MS2_infos_grip( GLMprm );
fprintf('\n \n');

n_prm_stroop_per_run = which_GLM_MS2_infos_stroop( GLMprm );
fprintf('\n \n');

switch gal.add_drv
    case 0
        drv_multiply = 1;
    case 1 % temporal derivative
        drv_multiply = 2;
    case 2 % temporal and spatial derivative
        drv_multiply = 3;
end

n_prm.per_run.RL        = n_prm_RL_per_run*drv_multiply;
n_prm.per_run.grip      = n_prm_grip_per_run*drv_multiply;
n_prm.per_run.stroop    = n_prm_stroop_per_run*drv_multiply;

% take into account movement regressors + presence or not of missed trials
% for global number of regressors

end % function