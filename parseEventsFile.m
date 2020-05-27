%% Parse Events File
%   Editor:    Daniel Elbich
%   Updated:   5/27/20
%
%   Script designed to parse an output containing event information, such
%   as onsets and timing, into separate files for each run. Files will be
%   saved as comma-separated values (CSV).
%
%   New behavioral file location will reflect orgnization below:
%
%   /StudyDir/rawdata/sub-s001/func/sub-s001_task-enc_run-1_events.csv
%   /StudyDir/rawdata/sub-s001/func/sub-s001_task-enc_run-2_events.csv
%   /StudyDir/rawdata/sub-s002/func/sub-s002_task-enc_run-1_events.csv
%   /StudyDir/rawdata/sub-s002/func/sub-s002_task-enc_run-2_events.csv
%
%   This is meant to be edited as needed on user end. Updates may occur to
%   the package at large but overall this only meant to standardized
%   required aspects of data (i.e. trial number, onsets, durations, conditions).

%% Set Analysis Parameters & Paths
% Load all relevent project information
if exist('flag','var') == 0
    
    %Select parameter file is flag does not exist
    [file,path]=uigetfile('*.mat','Select params file');
    filename=fullfile(path,file);
    load(filename);
    
end

%% Main Code

for i = 1:length(subjects)
    
    % Creates path to the current subjects behavioral file
    curSubj.behavDir  = [directory.Project filesep rawData.behavDir filesep...
        subjects{i} filesep 'behav'];
    curSubj.behavFile = dir([curSubj.behavDir filesep ...
        '*' taskInfo.Name '*.' rawData.behavFile]);
    
    if length(curSubj.behavFile)==1
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
    for curRun = 1:taskInfo.Runs
        
        %-- Initialize the names, onsets, durations, and pmods structure arrays
        % This section preallocates the names, onsets, and durations
        % structure arrays.
        
        % Convert raw onset column from msec to sec (divide by 1000)
        onsets = num2cell(BehavData.Onset(BehavData.Run == curRun)/1000)';
        
        % Set trial duration to zero for each trial, a stick function
        numTrials = length(onsets);
        durations = num2cell(zeros(1,numTrials));
        
        % Initialize cell of trial names
        currRunIDs = find(BehavData.Run == curRun);
        names = cell(1,numTrials);
        
        % Loop over all the trials in current run
        for ii = 1:numTrials
            
            if ii==1
                runEventData{1,1}='Trial';
                runEventData{1,2}='Condition';
                runEventData{1,3}='Onset';
                runEventData{1,4}='Duration';
                runEventData{1,5}='Run';  
            end
            
            runEventData{ii+1,1}=num2str(ii);
            runEventData{ii+1,2}=BehavData.EncodingCond{currRunIDs(ii)};
            runEventData{ii+1,3}=BehavData.Onset(currRunIDs(ii));
            runEventData{ii+1,4}=BehavData.ResponseTime(currRunIDs(ii));
            runEventData{ii+1,5}=num2str(curRun);
        end
        
        % Save individual event file for current run
        file = fopen([directory.Project filesep rawData.behavDir filesep ...
            subjects{i} filesep 'func' filesep subjects{i} '_task-'...
            taskInfo.Name '_run-' num2str(curRun) '_events.csv'], 'w');
        
        for a=1:size(runEventData,1)
            for b=1:size(runEventData,2)
                var = eval('runEventData{a,b}');
                try
                    fprintf(file, '%s', var);
                end
                fprintf(file, ',');
            end
            fprintf(file, '\n');
        end
        fclose(file);
        
    end
    
end

clc;
clear;
disp('All finished!');

