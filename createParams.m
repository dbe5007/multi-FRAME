%% Create Parameters
%   Editor:    Daniel Elbich
%   Created:   3/14/19
%
%   Creates parameters .mat file for use with MVPA scripts
%
%
%   Updated:   5/22/19
%
%   Reorganized script to handle variable preprocessing.
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


%% Set Pipeline Parameters

%  Set Path Variables
projectPath = '/path/to/project/parent/folder';
funcDir = inputdlg(['Enter Path to Functional Data Starting from: ' projectPath],...
    'File Path',[1 35],{'e.g. CPM/subj/behav'});
funcDir = [projectPath filesep funcDir];

%trialTag.Dir = [serverPath filesep '/path/to/trialtags.xlsx'];
trialTag.Dir = inputdlg(['Enter Path to Behavioral File Starting from: ' projectPath],...
    'File Path',[1 35],{'e.g. CPM/subj/behav'})

trialTag.ext = inputdlg('Enter Trial Tag File Extension:',...
    'File extension',[1 35],{'e.g. xlsx, csv'});

%  Functional Data Preprocessing Program
if exist('commandFlag','var')==0
    preprocPipeline = questdlg('Select Preprocessing Pipeline',...
        'Confirm Pipeline',...
        'SPM12','fMRIPrep','Cancel','fMRIPrep');
else
    preprocPipeline = 'fMRIPrep';
end

%  Obtain Data Parameter Information

try
numDat.Datapoints = inputdlg('How many volumes [Length of task]?',...
    'Volume Number',[1 35],{'150'});
numDat.Datapoints = str2double(num.Datapoints{:});

numDat.Runs = inputdlg('How many runs [# of functional runs]?',...
   'Number of Runs',[1 35],{'1'});
numDat.Runs = str2double(num.Runs{:});

analysisName = inputdlg('Enter Task Name: [as written in functional filename]',...
  'Task Name',[1 35],{'e.g. encoding, nback, rest'});

conditions = inputdlg('Enter Conditions of Interest: [separated by commas]',...
'Conditions',[1 35],{'e.g. face, object, place'});

%% Set Project-Specific Variables

Analysis.name            = ['SingleTrialModel' analysisName];
Analysis.behav.directory = [serverPath filesep trialTag.Dir];

% Please specify:
% -The units used (i.e., 'scans' or 'secs')
% -The TR or repition time
Model.units = 'secs';
Model.TR    = 2.5;
Model.TR  = inputdlg('Enter TR: [in seconds]',...
   'Repetition Time',[1 35],{'2'});
Model.TR  = str2double(Model.TR {:});

% Please specify if a mask is used during model estimation
Mask.on   = 0; %Default - no mask used
Mask.dir  = '/path/to/mask/directory';
Mask.name = 'name_of_mask.nii';

switch preprocPipeline
case 'SPM12'

Func.dir         = [serverPath filesep funcDir];
Func.wildcard    = '^warun.*\.nii'; % File
Mot.dir          = [serverPath filesep funcDir];
Func.motwildcard = '^rp_.*\.txt';

case 'fMRIPrep'
Func.dir         = [serverPath filesep funcDir];
Func.wildcard    = '*\analysisName.*\.nii.gz'; % File
Mot.dir          = [serverPath filesep funcDir];
Func.motwildcard = ['*\analysisName.*\.' trialTag.ext];

end

catch
    warning('Unable to set project information!.')
end

%% Project Paths

% General path setup
analysisPath = [serverPath filesep 'multivariate'];
studyPath   = [analysisPath filesep 'models/' Analysis.name];
setenv('studyPath',studyPath);
!mkdir -p $studyPath

        % Select pattern classification analysis. No searchlight is default
        if exist('commandFlag','var')==0
            classType = questdlg('Select Multivariate Analysis',...
                'Confirm Classification',...
                'MVPA','RSA','ERS','Cancel','MVPA');
        else
            classType = 'MVPA';
        end
        
        % Select analysis type. No searchlight is default
        if exist('commandFlag','var')==0
            analysisType = questdlg('Select Analysis Level for Multivariate Test:',...
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
                roiPath     = [projectPath 'ROIs/searchlight'];
                searchlight = 'Yes';
            otherwise
                roiPath     = '/path/to/common/mask/folder';
                searchlight = 'No';
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
                %if exist('commandFlag','var')==0
                    %trialAnalysis = questdlg('Test/Train on Runs or Trials',...
                        %'Confirm Trial/Run',...
                        %'Run','Trial','Cancel','Run');
                %else
                    trialAnalysis = 'Run';
                %end
                
                %if strcmpi(trialAnalysis,'Trial')==1
                    %leaveRuns = questdlg('Number of Trials to Test',...
                        %'Confirm Trial Tests',...
                        %'One','Two','Cancel','Run');
                %else
                    leaveRuns = NaN;
                %end
                
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
                %if exist('commandFlag','var')==0
                    %trialAnalysis = questdlg('Individual or Mean Conditions',...
                        %'Confirm Trial/Run',...
                        %'Individual','Mean','Cancel','Individual');
                %else
                    trialAnalysis = 'Individual';
                %end
                
                    case 'ERS'
                        [index,tf] = listdlg('Name','Possible Tasks',...
                            'PromptString','Ctrl+Click to select conditions:',...
                            'ListString',tasks,...
                            'ListSize',[280,300]);
                        tasks=tasks(index);
                        
                        save([projectPath 'vars/tasks.mat'],'tasks');
                
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
    subjDir=dir(funcDir);
    
    % Remove non-subject directories & possible files in directory
    for i=1:length(subjDir)
        subFlag(i) = isempty(strfind(subjDir(i).name,'.'));
    end
    subjDir = subjDir(subFlag);
    
    % Search 'study_path' directory to get list of subjects
    for i=1:length(subjDir)
        subjects{i,1}=subjDir(i).name;
    end
    
    !mkdir -p $analysisPath/vars
    save([analysisPath 'vars/subjects.mat'],'subjects');
    
    clear subjCount i subjDir;
catch
    warning('Unable to create subject list. Set to debug mode.');
end

%% Assign Conditions

try
    if exist('commandFlag','var')==0

        conds = [cellstr(conditions)];
        
        subconditionFlag=questdlg('Are there subconditions? (If unsure select No)',...
            'Confirm Subconditions','Yes','No','Cancel','No');
        
        save([analysisPath 'vars/conds.mat'],'conds');
        [index,tf] = listdlg('Name','Possible Conditions',...
            'PromptString','Ctrl+Click to select conditions:',...
            'ListString',conds,...
            'ListSize',[280,300]);

        subconds=subconds(index);
        save([analysisPath 'vars/subconds.mat'],'subconds');

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
switch analysisType
case 'Searchlight'
filename=strcat(analysisPath,'/params_',analysisName,'_',classType,'_',analysisType,'_',metric,'_',num2str(searchlightSize),'_',conds{1,1},'_',conds{1,2},'_Bootstrap_',bootstrap.flag,'.mat');
otherwise
filename=strcat(analysisPath,'params_',analysisName,'_',classType,'_',analysisType,'_',conds{1,1},'_',conds{1,2},'_subConds',num2str(exist(subconds,1)),'_Bootstrap_',bootstrap.flag,'.mat');
end

catch
    warning('Unable to set params filename. Set to debug mode.');
end

%% Save Params File
        switch analysisType
            case 'Searchlight'
                save(filename,'projectPath','funcDir','trialTag',...
                    'preprocPipeline','numDat','analysisName','conditions','Analysis',...
                    'Model','Mask','Func','Mot','analysisPath',...
                    'studyPath','metric','classType','analysisType','regressRT','roiPath','searchlight','bootstrap','trialAnalysis','leaveRuns','metric','searchlightSize',);
            otherwise
                save(filename,'analysisName','projectPath','studyPath',...
                'roi_path','subjects','conds','rois','analysisType',...
                'analysisName','classType','trialAnalysis','leaveRuns',...
                'bootstrap','metric','searchlightSize','regressRT');
end

clear;
clc;
disp('All finished!!');


%%%%%%% DELETE AFTER HERE %%%%%%%%


        case 'RSA'
            % Compute RSA with mean activation pattern (average all trials)
            % or individual trial pattern (all trials are separate)
            try
                %if exist('commandFlag','var')==0
                    %trialAnalysis = questdlg('Individual or Mean Conditions',...
                        %'Confirm Trial/Run',...
                        %'Individual','Mean','Cancel','Individual');
                %else
                    trialAnalysis = 'Individual';
                %end
                
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
    subjDir=dir(funcDir);
    
    % Remove non-subject directories & possible files in directory
    for i=1:length(subjDir)
        subFlag(i) = isempty(strfind(subjDir(i).name,'.'));
    end
    subjDir = subjDir(subFlag);
    
    % Search 'study_path' directory to get list of subjects
    for i=1:length(subjDir)
        subjects{i,1}=subjDir(i).name;
    end
    
    !mkdir -p $analysisPath/vars
    save([analysisPath 'vars/subjects.mat'],'subjects');
    
    clear subjCount i subjDir;
catch
    warning('Unable to create subject list. Set to debug mode.');
end

%% Assign Conditions

try
    if exist('commandFlag','var')==0

        conds = [cellstr(conditions)];
        
        subconditionFlag=questdlg('Are there subconditions? (If unsure select No)',...
            'Confirm Subconditions','Yes','No','Cancel','No');
        
        save([analysisPath 'vars/conds.mat'],'conds');
        [index,tf] = listdlg('Name','Possible Conditions',...
            'PromptString','Ctrl+Click to select conditions:',...
            'ListString',conds,...
            'ListSize',[280,300]);

        subconds=subconds(index);
        save([analysisPath 'vars/subconds.mat'],'subconds');

catch
    warning('Unable to assign conditions. Set to debug mode.');
end






