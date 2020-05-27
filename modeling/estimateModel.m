%% EstimateModel
%   Editor:    Daniel Elbich
%   Updated:   2/28/19
%
%   Script designed to for estimating a GLM specified using the
%   SpecifyModel.m script.  Allows user to display the trial type and
%   onsets/durations in the SPM Batch GUI.
%
% See also:  SpecifyModel, createParams

%% Set Analysis Parameters & Paths
% Load all relevent project information
if exist('flag','var') == 0
    
    %Select parameter file is flag does not exist
    [file,path]=uigetfile('*.mat','Select params file');
    filename=fullfile(path,file);
    load(filename);
    
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

spm('Defaults','FMRI')
spm_jobman('initcfg')
clc;

for curSub = 1:length(subjects) %for curSub = number %1:length(Subjects)
    
    fprintf('\n')
    fprintf('Subject: %s\n\n',subjects{curSub})
    
    % Model Directory: directory containing this subject's model
    Model.directory = fullfile(directory.Model, subjects{curSub});
    
    % If we are using a mask, create a path to the mask
    if Mask.on == 1
        Model.mask  = fullfile(Mask.dir, subjects{curSub}, Mask.name);
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
    
    % Get directory of motion files
    motionFiles = dir(fullfile(directory.Model, subjects{curSub},...
                ['*' taskInfo.Name '*.txt']));
    
    switch preprocPipeline
        case 'spm12'
            
            procDataDir = fullfile(directory.Project, 'derivatives',...
                'spmPreprocessing');
            
            procFuncFiles = dir(fullfile(procDataDir, subjects{curSub},...
                'func', 'run*',['ar*' taskInfo.Name '*.nii.gz']));
            
            for i = 1:taskInfo.Runs
                % Gunzip functional file
                
                setenv('modelData',[procFuncFiles(i).folder filesep procFuncFiles(i).name]);
                !gunzip $modelData
                
                %fprintf('Run: %d\n', i)
                %curFuncDir = fullfile(procDataDir, subjects{curSub}, ['func/run' num2str(i)]);
                %curMotDir = fullfile(Mot.dir, [Subjects{curSub} '_spm12'], [Runs{i} suffix]);
                
                %Model.runs{i}.scans = cellstr(spm_select('ExtFPList', curFuncDir, Func.wildcard, Inf));
                Model.runs{i}.scans = cellstr(spm_select('ExtFPList', ...
                    procFuncFiles(i).folder, Func.wildcard, Inf));
                Model.runs{i}.multicond = fullfile(Model.directory, SpecModelMats{i});        % from Model Spec
                Model.runs{i}.motion = [motionFiles(i).folder filesep motionFiles(i).name];
            end
            
            
        case 'fmriprep'
            
            procFuncFiles = dir(fullfile(directory.Project, 'preprocessing',...
                'fmriprep', subjects{curSub}, 'func', ...
                ['*' taskInfo.Name '*_bold.nii.gz']));
            
            for i = 1:taskInfo.Runs
                curFuncDir = fullfile(directory.Project, 'preprocessing',...
                    'fmriprep', subjects{curSub}, 'func');
                
                % Copy raw functional to Model Directory and gunzip
                setenv('data',[procFuncFiles(i).folder filesep procFuncFiles(i).name]);
                setenv('dest',[fullfile(directory.Model, subjects{curSub})]);
                setenv('modelData',[fullfile(directory.Model, subjects{curSub},...
                    procFuncFiles(i).name)]);
                !cp $data $dest
                !gunzip $modelData
            
                Model.runs{i}.scans = cellstr(spm_select('ExtFPList', ...
                    fullfile(directory.Model, subjects{curSub}), ...
                    [Func.wildcard num2str(i)], Inf));
                Model.runs{i}.multicond = fullfile(Model.directory, SpecModelMats{i});        % from Model Spec
                Model.runs{i}.motion = [motionFiles(i).folder filesep motionFiles(i).name]; % from realignment
            
            end
            
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
            fprintf('ERROR ON: %s', subjects{curSub})
        end
    end
    
    %% Manage resulting files
    switch preprocPipeline
        case 'spm12'
            % Set Gzip directories
            subjFuncDir = fullfile(directory.Model, subjects{curSub});
            setenv('procDataDir',[procDataDir filesep subjects{curSub}]);
            setenv('subjFuncDir',subjFuncDir);
            
            % Gzip output model
            !gzip $subjFuncDir/*.nii
            !gzip $procDataDir/func/*/*.nii
            
        case 'fmriprep'
            % Remove copied functional data from model directory
            !rm $dest/*_bold.nii
            
            % Gzip output model
            !gzip $dest/*.nii
            
    end
    
    % Create new SPM mat file for gzip nifti files
    load([directory.Model filesep subjects{curSub} filesep 'SPM.mat']);
    for i=1:length(SPM.Vbeta)
        SPM.Vbeta(i).fname = strrep(SPM.Vbeta(i).fname,'nii','nii.gz');
    end
    save([directory.Model filesep subjects{curSub} filesep 'SPM_gz.mat'],'SPM');
    
    clear SpecModelMats NumOfRuns curFuncDir curMotDir matlabbatch SPM;
    Model = rmfield(Model,'runs');
    Model = rmfield(Model,'directory');
    
    
end

clear;
clc;
disp('All finished!!');
