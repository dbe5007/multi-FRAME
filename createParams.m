%% Create params
%   Editor:    Daniel Elbich
%   Created:   3/14/19
%
%   Creates parameters .mat file for use with MVPA scripts
%
%
%   Updated:   9/19/19
%
%   Major updates:
%   -Added flag for regressing reaction time (RT) from betas prior to
%   classification
%   -Added flag and variables for bootstrap MVPA classification
%   -Added ERS functionality
%   
%
%   Minor updates:
%   -Cleaned up code to get subject list/directories
%   -Cleaned up code saving mat files with ROI names to preven overwriting
%   between ROI & Searchlight analyses
%
%
%   Updated:   4/22/19
%
%   Major updates:
%   -Added section to create separate params file for use with Specify
%   Model script. This will be deleted during creation of final param mat
%   file.
%   -Rescling now completed through AFNI. Must loaded/added to system PATH
%
%   Minor updates:
%   -Searches common folder for ROI masks instead of project folder
%   -List dialog to select 1 or more regions for reslicing
%   -User input to select searchlight/no searchlight anlaysis
%   -User input to select conditions of interest
%   -Saves subjects, rois, and conds variables to var directory in project
%   folder


%% Set Project-Specific Variables

try
    curProjects = {'FaceScene','ICEE','MORF','MAPP','FAME8','NEWPROJECTTEMPLATE'};
    %retiredProjects = {};
    
    index = listdlg('PromptString','Select Project:',...
        'SelectionMode','single',...
        'ListString',curProjects);
    
    switch curProjects{index}
        case 'FaceScene'
            if exist('commandFlag','var')==0
                analysisName = questdlg('Encoding or Retrieval','Task',...
                    'Encoding','Retrieval','Cancel','Encoding');
            else
                analysisName = 'Encoding';
            end
            
            tasks={'Encoding','Retrieval'};
            conds={'Target','Lure','Trained','Novel'};
            numRows = 195;
            serverPath = '/gpfs/group/nad12/default/nad12/facescene';
            
            switch analysisName
                case 'Encoding'
                    funcDir = 'Func_enc';
                    Analysis.behav.regexp = '/*ENCdm.wxls';
                    Number.OfRows = 195;
                    suffix = 'encoding';
                    runs = 5;
                case 'Retrieval'
                    funcDir = 'Func_ret_new';
                    Analysis.behav.regexp = '/*ret.xls';
                    Number.OfRows = 220;
                    suffix = 'ret';
                    runs = 5;
            end
            
            % Set Trial Duration in seconds
            regressRT.trialSec = 4;
            
            % Parameters for model estimation
            Analysis.name            = ['SingleTrialModel' analysisName];
            Analysis.behav.directory = [serverPath filesep 'Behav'];
            
            Func.dir         = [serverPath filesep funcDir];
            Func.wildcard    = '^warun.*\.nii'; % File
            Mot.dir          = [serverPath filesep funcDir];
            Func.motwildcard = '^rp_.*\.txt';
            
            % Please specify:
            % -The units used (i.e., 'scans' or 'secs')
            % -The TR or repition time
            Model.units = 'secs';
            Model.TR    = 2.5;
            
            % Please specify if a mask is used during model estimation
            Mask.on   = 0; %Default - no mask used
            Mask.dir  = '/path/to/mask/directory';
            Mask.name = 'name_of_mask.nii';
            
            % Get list of subject directories
            folders = dir([serverPath filesep funcDir]);
            for i=1:length(folders)
                subFlag(i) = ~isempty(strfind(folders(i).name,'_spm12'));
            end
            
            folders = folders(subFlag);
            for i=1:length(folders)
                Subjects{i,1}=erase(folders(i).name,'_spm12');
            end
            
        case 'ICEE'
            if exist('commandFlag','var')==0
                analysisName = questdlg('Encoding or Retrieval','Task',...
                    'Encoding','Retrieval','Cancel','Encoding');
            else
                analysisName = 'Encoding';
            end
            
            tasks={'Encoding','Retrieval'};
            conds={'cong', 'inc'};
            numRows = 195;
            serverPath = '/gpfs/group/nad12/default/nad12/ICEE';
            
            switch analysisName
                case 'Encoding'
                    funcDir = 'Func_enc';
                    Analysis.behav.regexp = '/*ENCdm.xls';
                    Number.OfRows = 195;
                    suffix = 'encoding';
                    runs = 5;
                case 'Retrieval'
                    funcDir = 'Func_ret_new';
                    Analysis.behav.regexp = '/*ret.xls';
                    Number.OfRows = 220;
                    suffix = 'ret';
                    runs = 5;
            end
            
            % Set Trial Duration in seconds
            regressRT.trialSec = 4;
            
            % Parameters for model estimation
            Analysis.name            = ['SingleTrialModel' analysisName];
            Analysis.behav.directory = [serverPath filesep 'Behav'];
            
            Func.dir         = [serverPath filesep funcDir];
            Func.wildcard    = '^warun.*\.nii'; % File
            Mot.dir          = [serverPath filesep funcDir];
            Func.motwildcard = '^rp_.*\.txt';
            
            % Please specify:
            % -The units used (i.e., 'scans' or 'secs')
            % -The TR or repition time
            Model.units = 'secs';
            Model.TR    = 2.5;
            
            % Please specify if a mask is used during model estimation
            Mask.on   = 0; %Default - no mask used
            Mask.dir  = '/path/to/mask/directory';
            Mask.name = 'name_of_mask.nii';
            
            % Get list of subject directories
            folders = dir([serverPath filesep funcDir]);
            for i=1:length(folders)
                subFlag(i) = ~isempty(strfind(folders(i).name,'_spm12'));
            end
            
            folders = folders(subFlag);
            for i=1:length(folders)
                Subjects{i,1}=erase(folders(i).name,'_spm12');
            end
            
        case 'MAPP'
            
            
    end
catch
    warning('Unable to set project information!.')
end


%% Project Paths

% General path setup
projectPath = [serverPath filesep 'RSA'];
studyPath   = [projectPath filesep 'models/unsmoothed/SingleTrialModel' analysisName];
setenv('project_path',projectPath);
!mkdir -p $project_path

% Select analysis type. No searchlight is default
if exist('commandFlag','var')==0
    stepFlag = questdlg('Has Specify/Estimate Model been run?','Confirm Step',...
        'Yes','No','Cancel','Yes');
else
    stepFlag = 'Yes';
end

switch stepFlag
    case 'Yes'
        
        % Select pattern classification analysis. No searchlight is default
        if exist('commandFlag','var')==0
            classType = questdlg('Select Classification Type',...
                'Confirm Classification',...
                'MVPA','RSA','Cancel','MVPA');
        else
            classType = 'MVPA';
        end
        
        % Select analysis type. No searchlight is default
        if exist('commandFlag','var')==0
            analysisType = questdlg('Select Analysis Type',...
                'Confirm Searchlight',...
                'Searchlight','ROI','Cancel','ROI');
        else
            analysisType = 'ROI';
        end
        
        % Account for RT/regress out. No is default
        if exist('commandFlag','var')==0
            regressRT.flag = questdlg('Regress out Reaction Time (RT)?',...
                'Confirm Regression',...
                'Yes','No','Cancel','No');
            
        else
            regressRT.flag = 'No';
        end
        
        switch analysisType
            case 'Searchlight'
                roi_path    = [projectPath 'ROIs/searchlight'];
                searchlight = 'Yes';
            otherwise
                roi_path    = '/path/to/common/mask/folder';
                searchlight = 'No';
        end
end

%% Params for Specify/Estimate Model

try
    switch stepFlag
        case 'No'
            
            for i=1:runs
                Runs{i} = ['run' num2str(i)];
            end
            
            % Parameters for model estimation
            Analysis.directory = studyPath;
            
            filename=strcat(projectPath,'specify_',analysisName,'_model_params2.mat');
            save(filename,'analysisName','projectPath','studyPath',...
                'serverPath','Analysis','funcDir','Subjects','Func',...
                'Mot','suffix','Runs','Model','Mask');
            
            if ~exist('commandFlag','var')
                msgbox(['Params for Specify Model script have been created. Please '...
                    'run Specify Model and then re-run this script.']);
            end
            
            clear;
            return;
        case 'Yes'
            %delete([project_path,'specify_',analysisName,'_model_params.mat']);
    end
catch
    warning(['Unable to create parameter file for Specify Model.'...
        'Set to debug mode.']);
end

%% Classification Flags

try
    % Bootstrapping Setup
    
    if exist('commandFlag','var')==0
        bootstrap.flag  = questdlg('Perform Bootstrap',...
            'Confirm Bootstrap',...
            'Yes','No','Cancel','No');
        
        if strcmpi(bootstrap.flag,'Yes')==1
            bootstrap.numRuns     = runs;
            bootstrap.numRows     = numRows;
            bootstrap.numTrials   = bootstrap.numRows/...
                bootstrap.numRuns;
            bootstrap.perm        = inputdlg('How many permutations?',...
                'Permutation Number',[1 35],{'1000'});
            bootstrap.perm = str2double(searchlightSize{:});
            
            bootstrap.trialsPerRun = randperm(length...
                (1:bootstrap.numTrials));
            
            for i=1:bootstrap.numRuns
                bootstrap.structNames{i} = ['perm' num2str(i)];
            end
        end
    end
    
    switch classType
        % MVPA Flags
        case 'MVPA'
            % Compute MVPA with entire run to train/test or individual
            % trials to train/test
            try
                if exist('commandFlag','var')==0
                    trialAnalysis = questdlg('Test/Train on Runs or Trials',...
                        'Confirm Trial/Run',...
                        'Run','Trial','Cancel','Run');
                else
                    trialAnalysis = 'Run';
                end
                
                if strcmpi(trialAnalysis,'Trial')==1
                    leaveRuns = questdlg('Number of Trials to Test',...
                        'Confirm Trial Tests',...
                        'One','Two','Cancel','Run');
                else
                    leaveRuns = NaN;
                end
                
                if strcmpi(analysisType, 'Searchlight')==1
                    metric = questdlg('Searchlight Metric',...
                        'SearchlightSize','count','radius','Cancel','count');
                    
                    searchlightSize = inputdlg('Size of the searchlight',...
                        'Searchlight Size',[1 35],{'50'});
                    searchlightSize = str2double(searchlightSize{:});
                end
                
            catch
                warning(['Error in setting MVPA analysis flags. '...
                    'Set to debug mode.']);
            end
            
            % RSA Flags
        case 'RSA'
            % Compute RSA with mean activation pattern (average all trials)
            % or individual trial pattern (all trials are separate)
            try
                if exist('commandFlag','var')==0
                    trialAnalysis = questdlg('Individual or Mean Conditions',...
                        'Confirm Trial/Run',...
                        'Individual','Mean','Cancel','Individual');
                else
                    trialAnalysis = 'Individual';
                end
                
                % Perform RSA or ERS. Default is RSA.
                if exist('commandFlag','var')==0
                    classificationType = questdlg...
                        ('Single Task RSA or ERS?','Confirm Similarity',...
                        'RSA','ERS','Cancel','RSA');
                else
                    classificationType = 'RSA';
                end
                
                switch classificationType
                    case 'ERS'
                        [index,tf] = listdlg('Name','Possible Tasks',...
                            'PromptString','Ctrl+Click to select conditions:',...
                            'ListString',tasks,...
                            'ListSize',[280,300]);
                        tasks=tasks(index);
                        
                        save([projectPath 'vars/tasks.mat'],'tasks');
                end
                
            catch
                warning(['Error in setting RSA/ERS analysis flags. '...
                    'Set to debug mode.']);
            end
    end
catch
    warning('Error in setting classification flags. Set to debug mode.');
end

%% Create Subject List

try
    % List subject directory
    subjDir=dir(studyPath);
    
    % Remove non-subject directories & possible files in directory
    for i=1:length(subjDir)
        subFlag(i) = isempty(strfind(subjDir(i).name,'.'));
    end
    subjDir = subjDir(subFlag);
    
    % Search 'study_path' directory to get list of subjects
    for i=1:length(subjDir)
        subjects{i,1}=subjDir(i).name;
    end
    
    !mkdir -p $project_path/vars
    save([projectPath 'vars/subjects.mat'],'subjects');
    
    clear subjCount i subjDir;
catch
    warning('Unable to create subject list. Set to debug mode.');
end

%% Assign Conditions

try
    if exist('commandFlag','var')==0
        conditionFlag=questdlg('Conditions specified?','Confirm Conditions',...
            'Yes','No','Cancel','Yes');
        
        subconditionFlag=questdlg('Are there subconditions? (If unsure select No)',...
            'Confirm Subconditions','Yes','No','Cancel','No');
        
        switch conditionFlag
            case 'Yes'
                load([projectPath 'vars/conds.mat']);
            case 'No'
                [index,tf] = listdlg('Name','Possible Conditions',...
                    'PromptString','Ctrl+Click to select conditions:',...
                    'ListString',conds,...
                    'ListSize',[280,300]);
                conds=conds(index);
                
                save([projectPath 'vars/conds.mat'],'conds');
        end
        
        switch subconditionFlag
            case 'Yes'
                load([projectPath 'vars/subConds.mat']);
        end
        
    else
        load([projectPath 'vars/conds.mat']);
        if exist([projectPath 'vars/subConds.mat'],'file')
            load([projectPath 'vars/subConds.mat']);
        end
    end
catch
    warning('Unable to assign conditions. Set to debug mode.');
end

%% Create ROI List

try
    % Flag to reslice regions/load in previously resliced regions. No is
    % default
    
    % Interactive run asks user input - command line defaults to reslice if roi
    % mat file is missing
    if ~exist('commandFlag','var')
        resliceFlag=questdlg('Are ROIs in subject space?','Confirm Reslicing',...
            'Yes','No','Cancel','No');
    else
        switch searchlight
            case 'Yes'
                if ~exist([projectPath 'vars/rois_searchlight.mat'],'file')
                    resliceFlag='No';
                else
                    resliceFlag='Yes';
                end
            otherwise
                if ~exist([projectPath 'vars/rois.mat'],'file')
                    resliceFlag='No';
                else
                    resliceFlag='Yes';
                end
        end
    end
    
    switch resliceFlag
        case 'Yes'
            switch searchlight
                case 'Yes'
                    load([projectPath 'vars/rois_searchlight.mat']);
                otherwise
                    load([projectPath 'vars/rois.mat']);
            end
        case 'No'
            
            % Get subject directories
            % dataDir=[study_path filesep subjects{1} filesep 'beta_0001.nii'];
            dataDir = [serverPath funcDir filesep subjects{1}...
                '_spm12' filesep 'run1encoding' filesep 'warun1encoding.nii'];
            
            % Load in lab ROIs
            maskFolders=dir(roi_path);
            maskFolders(1:2)=[];
            maskCount=1;
            
            % Get lists of possible masks in common LabMasks folder
            switch analysisType
                case 'Searchlight'
                    % Set prefix for save path
                    roiPrefix='ROIs/searchlight/reslice_';
                    
                    for i=1:length(maskFolders)
                        if maskFolders(i).isdir==0
                            maskName{maskCount}=maskFolders(i).name;
                            maskList{maskCount}=[maskFolders(i).folder...
                                filesep maskFolders(i).name];
                            maskCount=maskCount+1;
                        end
                    end
                otherwise
                    % Set prefix for save path
                    roiPrefix='ROIs/reslice_';
                    
                    % Search through all folders and create list of masks
                    for i=1:length(maskFolders)
                        
                        folder=dir([maskFolders(i).folder filesep maskFolders(i).name]);
                        folder(1:2)=[];
                        
                        for ii=1:length(folder)
                            if folder(ii).isdir==0
                                maskName{maskCount}=folder(ii).name;
                                maskList{maskCount}=[folder(ii).folder...
                                    filesep folder(ii).name];
                                maskCount=maskCount+1;
                            end
                        end
                    end
            end
            
            % Drop-down list of found masks - select all that apply
            [index,tf] = listdlg('Name','Available Masks',...
                'PromptString','Ctrl+Click to choose masks:',...
                'ListString',maskName,...
                'ListSize',[280,300]);
            maskList=maskList(index);
            regions=maskName(index);
            
            for i=1:length(regions)
                datainput=maskList{i};
                reference=dataDir;
                output=strcat(projectPath,roiPrefix,regions{i});
                
                % Set system/terminal variables
                setenv('region',regions{i});
                setenv('datainput',datainput);
                setenv('reference',reference);
                setenv('output',output);
                
                % Skip previously resliced regions - AFNI will NOT
                % overwrite!!!
                if ~exist(output,'file')
                    % Call AFNI from system/terminal to fit region to subject space
                    !echo "Fitting $region to subject functional..."
                    !3dresample -master $reference -prefix $output -input $datainput
                end
                
                rois{1,i}=strcat('reslice_',regions{i});
            end
            
    end
    
    roi_path = [projectPath 'ROIs'];
    
    % Save roi list to separate .mat file
    switch searchlight
        case 'Yes'
            save([projectPath 'vars/rois_searchlight.mat'],'rois');
        otherwise
            save([projectPath 'vars/rois.mat'],'rois');
    end
catch
    warning('Unable to create region list. Set to debug mode.');
end

%% Set filename for output

try
    % Set filename for parameter .mat file
    switch classType
        case 'MVPA'
            switch analysisType
                case 'Searchlight'
                    filename=strcat(projectPath,'params_',...
                            analysisName,'_',classType,'_',...
                            analysisType,'_',trialAnalysis,'_',...
                            metric,'_',num2str(searchlightSize),'_',...
                            conds{1,1},'_',conds{1,2},'_Bootstrap_',...
                            bootstrap.flag,'.mat');
                otherwise
                    if exist('subConds','var')
                        filename=strcat(projectPath,'params_',...
                            analysisName,'_',classType,'_',...
                            analysisType,'_',trialAnalysis,'_',...
                            conds{1,1},'_',conds{1,2},'_',...
                            subConds{1,1},'_',subConds{1,2},'_Bootstrap_',...
                            bootstrap.flag,'.mat');
                    else
                        filename=strcat(projectPath,'params_',...
                            analysisName,'_',classType,'_',...
                            analysisType,'_',trialAnalysis,'_',...
                            conds{1,1},'_',conds{1,2},'_Bootstrap_',...
                            bootstrap.flag,'.mat');
                    end
            end
        case 'RSA'
            if exist('subConds','var')
                filename=strcat(projectPath,'params_',...
                    analysisName,'_',classificationType,'_',...
                    analysisType,'_',trialAnalysis,'_',...
                    conds{1,1},'_',conds{1,2},'_',...
                    subConds{1,1},'_',subConds{1,2},'_Bootstrap_',...
                    bootstrap.flag,'.mat');
            else
                
                filename=strcat(projectPath,'params_',...
                    analysisName,'_',classificationType,'_',...
                    analysisType,'_',trialAnalysis,'_',...
                    conds{1,1},'_',conds{1,2},'_Bootstrap_',...
                    bootstrap.flag,'.mat');
            end
    end
catch
    warning('Unable to set params filename. Set to debug mode.');
end

%% Save Params File
switch classType
    case 'MVPA'
        switch analysisType
            case 'Searchlight'
                save(filename,'analysisName','projectPath','studyPath',...
                    'roi_path','subjects','conds','rois','analysisType',...
                    'analysisName','classType','trialAnalysis','leaveRuns',...
                    'bootstrap','metric','searchlightSize','regressRT');
            otherwise
                if exist('subConds','var')
                    save(filename,'analysisName','projectPath','studyPath',...
                        'roi_path','subjects','conds','subConds','rois',...
                        'analysisType','analysisName','classType',...
                        'bootstrap','trialAnalysis','leaveRuns','regressRT');
                else
                    save(filename,'analysisName','projectPath','studyPath',...
                        'roi_path','subjects','conds','rois','analysisType',...
                        'analysisName','classType','trialAnalysis',...
                        'bootstrap','leaveRuns','regressRT');
                end
        end
    case 'RSA'
        if exist('subConds','var')
            save(filename,'analysisName','projectPath','studyPath',...
                'roi_path','subjects','conds','subConds','rois',...
                'analysisType','analysisName','classType','bootstrap',...
                'trialAnalysis','classificationType','regressRT');
        else
            save(filename,'analysisName','projectPath','studyPath',...
                'roi_path','subjects','conds','rois','analysisType',...
                'analysisName','classType','trialAnalysis','bootstrap',...
                'classificationType','tasks','regressRT');
        end
end

clear;
clc;
disp('All finished!!');
