%% EstimateModel
%   Editor:    Daniel Elbich
%   Updated:   2/28/19
%
%   Script designed to for estimating a GLM specified using the
%   SpecifyModel.m script.  Allows user to display the trial type and
%   onsets/durations in the SPM Batch GUI.
%
% See also:  SpecifyModel, createParams

%% User Input
% You should ONLY (!!!!!!) need to edit this highlighted section of the
% script.

% Add SPM12 to path
addpath(genpath('/path/to/spm12folder'));

% User Input Step 1: Analysis, Funcs, and Masks

% Set Analysis Parameters & Paths
% Load subject IDs, ROIs, and Condition flags
if exist('flag','var') == 0
    
    %Select parameter file is flag does not exist
    uiopen('*.mat')
    
end

% User Input Step 1: Options
% Set the following jobman_option to 'interactive' to view in SPM parameters the GUI.
% Press any key into the command window to continue to next one sample t test.
% Set the following jobman option to 'run' to skip the viewing of the
% SPM parameters in the GUI and go directly to running of the one
% sample t-tests

show          = 0; % 0 = input as multiple conditions file, 1 = input parameters for showing in GUI
jobman_option = 'run'; % interactive = show in GUI, run = run through SPM

%% Routine

clc
fprintf('Analysis: %s\n\n', Analysis.name)
fprintf('Gathering Data...\n\n')
spm('Defaults','FMRI')
spm_jobman('initcfg')

for curSub = 1:length(Subjects) %for curSub = number %1:length(Subjects)
    
    fprintf('\n')
    fprintf('Subject: %s\n\n',Subjects{curSub})
    
    % Model Directory: directory containing this subject's model
    Model.directory = fullfile(Analysis.directory, Subjects{curSub});
    
    % If we are using a mask, create a path to the mask
    if Mask.on == 1
        Model.mask  = fullfile(Mask.dir, Subjects{curSub}, Mask.name);
    end
    
    % Find the SpecModel *.mat files. These should be in the model
    % directory, defined above
    SpecModelMats   = cellstr(spm_select('List', Model.directory, '.*Run.*\.mat'));
    
    % If we do not find any *.mat files, give an error informing the user
    % that this has occured
    if cellfun('isempty', SpecModelMats)
        error('Could not find Specify Model *.mat files') %#ok<*NODEF>
    end
    
    % Determine number of runs from number of Model Spec *.mat files found
    NumOfRuns = length(SpecModelMats);
    
    for i = 1:NumOfRuns
        fprintf('Run: %d\n', i)
        curFuncDir              = fullfile...
            (Func.dir, [Subjects{curSub} '_spm12'], [Runs{i} suffix]);
        curMotDir               = fullfile...
            (Mot.dir, [Subjects{curSub} '_spm12'], [Runs{i} suffix]);
        Model.runs{i}.scans     = cellstr...
            (spm_select('ExtFPList', curFuncDir, Func.wildcard, Inf));
        Model.runs{i}.multicond = fullfile...
            (Model.directory, SpecModelMats{i});        % from Model Spec
        Model.runs{i}.motion    = spm_select...
            ('FPList', curMotDir, Func.motwildcard); % from realignment
    end
    
    %% Set and Save the SPM job
    
    try
    onsets    = [];
    durations = [];
    names     = [];
    pmod      = [];
    
    % Directory
    matlabbatch{1}.spm.stats.fmri_spec.dir = {Model.directory};
    
    % Model Parameters
    matlabbatch{1}.spm.stats.fmri_spec.timing.units   = Model.units;
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT      = Model.TR;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t  = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
    
    % Session Specific
    for curRun = 1:length(Model.runs)
        matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).scans = Model.runs{curRun}.scans; %#ok<*AGROW>
        if show == 0
            matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).cond  = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).multi = {Model.runs{curRun}.multicond};
        elseif show == 1
            load(Model.runs{curRun}.multicond)
            for curTrialType = 1:length(names)
                matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).cond(curTrialType).name     = names{curTrialType};
                matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).cond(curTrialType).onset    = onsets{curTrialType};
                matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).cond(curTrialType).duration = durations{curTrialType};
                matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).cond(curTrialType).tmod     = 0;
                if isempty(pmod)
                    matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).cond(curTrialType).pmod = struct('name', {}, 'param', {}, 'poly', {});
                else
                    matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).cond(curTrialType).pmod.name  = pmod(curTrialType).name{1};
                    matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).cond(curTrialType).pmod.param = pmod(curTrialType).param{1};
                    matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).cond(curTrialType).pmod.poly  = pmod(curTrialType).poly{1};
                end
            end
        end
        matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).regress   = struct('name', {}, 'val', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).multi_reg = {Model.runs{curRun}.motion};
        matlabbatch{1}.spm.stats.fmri_spec.sess(curRun).hpf       = 128;
    end
    
    % Misc
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    if Mask.on == 1
        matlabbatch{1}.spm.stats.fmri_spec.mask = {Model.mask};
    else
        matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    end
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    % Estimation
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep;
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tname = 'Select SPM.mat';
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).name  = 'filter';
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).value = 'mat';
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).name  = 'strtype';
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1).sname = 'fMRI model specification: SPM.mat File';
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_output = substruct('.','spmmat');
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    catch
        disp('Unable to create matlabbatch!');
        fprintf('ERROR ON: %s', Subjects{curSub});
    end
    
    save(fullfile(Model.directory, 'Job.mat'), 'matlabbatch');
    clear onsets durations names pmod;
    
    %% Run the SPM job
    
    if strcmp(jobman_option,'interactive')
        fprintf('\n')
        fprintf('Displaying SPM Job...\n')
        spm_jobman(jobman_option, matlabbatch)
        pause
    elseif strcmp(jobman_option,'run')
        try
            spm_jobman(jobman_option, matlabbatch)
        catch ER %#ok<*NASGU>
            disp(ER)
            fprintf('ERROR ON: %s', Subjects{curSub})
        end
    end
    
    % Gzip all functionals
    subjFuncDir = fullfile(Func.dir, Subjects{curSub},'func');
    setenv('subjFuncDir',subjFuncDir);
    !gzip $subjFuncDir/*/*.nii
    
    % Create new SPM mat file for gzip nifti files
    load([studyPath filesep Subjects{curSub} filesep 'SPM.mat']);
    for i=1:length(SPM.Vbeta)
        SPM.Vbeta(i).fname = strrep(SPM.Vbeta(i).fname,'nii','nii.gz');
    end
    save([studyPath filesep Subjects{curSub} filesep 'SPM_gz.mat'],'SPM');
    
    clear SpecModelMats NumOfRuns curFuncDir curMotDir matlabbatch SPM;
    Model = rmfield(Model,'runs');
    Model = rmfield(Model,'directory');
    
    
end

disp('All finished!!');
