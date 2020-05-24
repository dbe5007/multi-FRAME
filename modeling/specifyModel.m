%% SpecifyModel
%   Editor:    Daniel Elbich
%   Updated:   2/27/19
%
%   Script designed to build multiple conditions files for later use with
%   SPM's matlabbatch system. Edited to inferface with project 'FaceScene'
%   file structure on Linux ACI system.
%
%   Assumes that the behavioral data are organized as follows:
%
%   /StudyDir/BehavDir/s001/ConcatenatedBehavioralData.csv
%   /StudyDir/BehavDir/s002/ConcatenatedBehavioralData.csv
%   /StudyDir/BehavDir/s003/ConcatenatedBehavioralData.csv
%
%   Where ConcatenatedBehavioralData.csv is a comma seperated text file
%   where all functional runs are concatenated, with a column that
%   indicates which rows (i.e., trials) belong to which functional run.
%
%   See also createParams.m, EstimateModel.m
%
%
%   Updates:
%   
%   5/6/19 - Number of trials, durations, and onsets saved as csv file for
%   data QA.

%% Parameter Setup

%==========================================================================
%                           User Input
%==========================================================================

% Set Analysis Parameters & Paths
% Load subjects, paths, and condition flags
if exist('flag','var')==0
    
    %Select parameter file interactively when running from MATLAB
    uiopen('*.mat')
    
end

%% Main Code
for i = 1:length(subjects)
    
    % Creates path to the current subjects behavioral file
    curSubj.behavDir  = [directory.Project filesep rawData.behavDir filesep...
        subjects{i} filesep 'behav'];
    curSubj.behavFile = dir([curSubj.behavDir filesep ...
        '*' taskName '*.' rawData.behavFile]);
    
    if length(length(curSubj.behavFile))==1
        % Reads in the subjects behavioral data using the readtable command.
        % See readtable for more details.
        fprintf('Reading in subject %s behavioral data ...\n', subjects{i});
        BehavData = readtable([curSubj.behavDir filesep curSubj.behavFile.name]);
        
        % Clean up variable names
        BehavData.Properties.VariableNames = regexprep(regexprep...
            (BehavData.Properties.VariableNames, '_', ''), '^x', '');
    end
    
    % Creates a path to this subjects analysis directory & creates that
    % directory if it does not already exist.
    curSubj.directory = fullfile(directory.Model, subjects{i});
    if ~isdir(curSubj.directory)
        mkdir(curSubj.directory)
    end
    
    % Initalize the counter cell array to track number of trials in each functional run
    %number_of_runs = max(unique(BehavData.block));
    fprintf('Sorting Behavioral Data...\n\n')
    
    % Build the multiple conditions *.mat file for each run
    for curRun = 1:dataInfo.Runs
        
        if length(length(curSubj.behavFile))>1
            % Reads in the subjects behavioral data using the readtable command.
            % See readtable for more details.
            fprintf('Reading in subject %s behavioral data ...\n', subjects{i});
            BehavData = readtable([curSubj.behavDir filesep curSubj.behavFile.name]);
            
            % Clean up variable names
            BehavData.Properties.VariableNames = regexprep(regexprep...
                (BehavData.Properties.VariableNames, '_', ''), '^x', '');
        end
        
        %-- Initialize the names, onsets, durations, and pmods structure arrays
        % This section preallocates the names, onsets, and durations
        % structure arrays.
        
        % Convert raw onset column from msec to sec (divide by 1000)
        onsets = num2cell(BehavData.Onset(BehavData.Run == curRun)/1000)';
        
        % Set trial duration to zero for each trial, a stick function
        number_of_trials = length(onsets);  % CAN DELETE ONE TRIAL NUMBER FROM CREATE PARAMS
        durations        = num2cell(zeros(1,number_of_trials));
        %durations        = num2cell(zeros(1,dataInfo.Trials));
        
        % Initialize cell of trial names
        currRunIDs = find(BehavData.Run == curRun);
        names   = cell(1,number_of_trials); % CAN DELETE ONE TRIAL NUMBER FROM CREATE PARAMS
        %names   = cell(1,dataInfo.Trials);
        
        % Loop over all the trials in current run
        for ii = 1:number_of_trials % CAN DELETE ONE TRIAL NUMBER FROM CREATE PARAMS
            % for ii = 1:dataInfo.Trials
            
            names{ii} = sprintf('trial-%s_condition-%s_retType-%s_dmScore-%s',...
                char(num2str(currRunIDs(ii))),...
                char(BehavData.EncodingCond(currRunIDs(ii))),...
                'Hit', char(BehavData.DMscore(currRunIDs(ii))));
            
        end
        
        %-- Save the Multiple Conditions *.mat file
        % Save the names, onsets, and durations variables to .mat file
        % to be used for later model estimation in SPM. See EstimateModel.m
        
        matfilename = fullfile(curSubj.directory, ...
            ['Run', num2str(curRun, '%03d'), '_multiple_conditions.mat']);
        fprintf('Saving subject %s run %d multiple conditions file...\n',...
            subjects{i}, curRun)
        save(matfilename, 'names', 'onsets', 'durations');
        
        % Summary Trial Numbers
        if i==1
            labelHeader={['names_' num2str(curRun, '%03d')],...
                ['onsets_' num2str(curRun, '%03d')],...
                ['durations_' num2str(curRun, '%03d')]};
            
            finalSummary{1,1}='Subject_ID';
            index=size(finalSummary,2);
            index=[index+1 index+3];
            finalSummary(1,index(1):index(2))=labelHeader;
        end
        
        if curRun==1
            finalSummary{i+1,1}=subjects{i};
            finalIndex=[2:4];
        end
        
        finalSummary{i+1,finalIndex(1)}=num2str(length(names));
        finalSummary{i+1,finalIndex(2)}=num2str(length(onsets));
        finalSummary{i+1,finalIndex(3)}=num2str(length(durations));
        
        finalIndex=finalIndex+3;
        
    end
    
    %ersTagFilename = fullfile(curSubj.directory,'ersTags.mat');
    %save(ersTagFilename,'ersTags');
        
end


%% Write output summary file

file = fopen([directory.Model filesep 'allTrialsSummary.csv'], 'w');

for a=1:size(finalSummary,1)
    for b=1:size(finalSummary,2)
        var = eval('finalSummary{a,b}');
        try
            fprintf(file, '%s', var);
        end
        fprintf(file, ',');
    end
    fprintf(file, '\n');
end

fclose(file);
clc;
disp('All Finished!!');
clear;
