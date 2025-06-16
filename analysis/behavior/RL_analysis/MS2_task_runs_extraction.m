function[runs] = MS2_task_runs_extraction(taskName, subject_id)
% [runs] = MS2_task_runs_extraction(taskName, subject_id)
% MS2_task_runs_extraction outputs the runs corresponding to the subject
% asked in the inputs (subject_id) for the corresponding task asked in the
% inputs (taskName)
% 
% INPUTS
% taskName: string containing task name
% 'grip', 'stroop' or 'RL'
%
% subject_id: subject identification name (string)
%
% OUTPUTS
% runs: corresponding runs for the subject*task asked in the inputs

switch subject_id
    
    % For the subjects who did the tasks in the following order: 1) RL, 2)
    % Stroop, 3) Grip, 4) RL, 5) Stroop, 6) Grip, 7) RL
    case {'sX_YYYYYY'} % replace with your own subject identifiers
        switch taskName
            case 'grip'
                runs = [3,6];
            case 'stroop'
                runs = [2,5];
            case 'RL'
                runs = [1,4,7];
        end
        
        % For the subjects who did the tasks in the following order: 1) RL, 2)
    % Grip, 3) Stroop, 4) RL, 5) Grip, 6) Stroop, 7) RL
    case {'sX_ZZZZZZ'} % replace with your own subject identifiers
        switch taskName
            case 'grip'
                runs = [2,5];
            case 'stroop'
                runs = [3,6];
            case 'RL'
                runs = [1,4,7];
        end
end

end