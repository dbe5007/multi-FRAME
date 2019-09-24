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

% Clean up and print update to the command window
clc;
fprintf('Model: %s\n\n', Analysis.name)
fprintf('Model Directory: \n')
disp(Analysis.directory)
fprintf('\n')
fprintf('Behavioral Data Directory: \n')
disp(Analysis.behav.directory)
fprintf('\n')

% Subject Loop
for i = 1:length(Subjects)
    
    % Creates path to the current subjects behavioral file
    curSubj.behavDir  = [Analysis.behav.directory '/' Subjects{i}];
    curSubj.behavFile = dir([curSubj.behavDir Analysis.behav.regexp]);
    
    % Reads in the subjects behavioral data using the readtable command.
    % See readtable for more details.
    fprintf('Reading in subject %s behavioral data ...\n', Subjects{i});
    BehavData = readtable([curSubj.behavDir filesep curSubj.behavFile.name]);
    
    % Clean up variable names
    BehavData.Properties.VariableNames = regexprep(regexprep...
        (BehavData.Properties.VariableNames, '_', ''), '^x', '');
    
    % Creates a path to this subjects analysis directory & creates that
    % directory if it does not already exist.
    curSubj.directory = fullfile(Analysis.directory, Subjects{i});
    if ~isdir(curSubj.directory)
        mkdir(curSubj.directory)
    end
    
    % Initalize the counter cell array to track number of trials in each functional run
    number_of_runs = max(unique(BehavData.block));
    fprintf('Sorting Behavioral Data...\n\n')
    
    % Build the multiple conditions *.mat file for each run
    for curRun = 1:number_of_runs
        
        %-- Initialize the names, onsets, durations, and pmods structure arrays
        % This section preallocates the names, onsets, and durations
        % structure arrays.
        
        % Convert raw onset column from msec to sec (divide by 1000)
        onsets = num2cell(BehavData.rawonset(BehavData.block == curRun)/1000)';
        
        % Set trial duration to zero for each trial, a stick function
        number_of_trials = length(onsets);
        durations        = num2cell(zeros(1,number_of_trials));
        
        % Initialize cell of trial names
        currRunIDs = find(BehavData.block == curRun);
        names   = cell(1,number_of_trials);
        
        % Loop over all the trials in current run
        for ii = 1:number_of_trials
            
            switch analysisName
                case 'TaskA'
                    try
                        %%% pull variables to add to the trial file name
                        % image names
                        
                        % change pic1 to pic for retrieval
                        imageFace = erase(BehavData.pic1(currRunIDs(ii)),'.\faces\');
                        imageScene = erase(BehavData.pic2(currRunIDs(ii)),'.\scenes\scenes\');
                        imageFace = erase(imageFace,'"');
                        imageScene = erase(imageScene,'"');
                        imageFace = strrep(imageFace,'''','');
                        imageScene = strrep(imageScene,'''','');
                        
                        % Merge image names
                        imagename = strcat(imageFace{1,1},'|',imageScene{1,1});
                        
                        % Study trial type
                        if BehavData.StudyTT(currRunIDs(ii)) == 0
                            studyType       = 'Trained';
                        elseif BehavData.StudyTT(currRunIDs(ii)) == 1
                            studyType       = 'Novel';
                        end
                        
                        % Retrieval trial type
                        if BehavData.RetTT(currRunIDs(ii)) == 0
                            retType         = 'Target';
                        elseif BehavData.RetTT(currRunIDs(ii)) == 1
                            retType         = 'Lure';
                        end
                        
                        % Difference-Memory score (DMscore)
                        switch BehavData.DMscore(currRunIDs(ii))
                            case 1 % Miss or False Alarm
                                if BehavData.RetTT(currRunIDs(ii)) == 0
                                    dmScore = 'Miss';
                                else
                                    dmScore = 'FalseAlarm';
                                end
                            case 2 % Familiar or Slight False Alarm
                                if BehavData.RetTT(currRunIDs(ii)) == 0
                                    dmScore = 'Familiar';
                                else
                                    dmScore = 'SlightFalseAlarm';
                                end
                            case 3 % Hit or Crrect Rejection
                                if BehavData.RetTT(currRunIDs(ii)) == 0
                                    dmScore = 'Hit';
                                else
                                    dmScore = 'CorrectRejection';
                                end
                            otherwise
                                dmScore = 'NoBehaviorRecorded';
                        end
                        
                        % informative, unique, BIDS style trial name
                        names{ii} = sprintf...
                            ('imagename-%s_studyType-%s_retType-%s_dmScore%s',...
                            imagename, studyType, retType, dmScore);
                        
                        ersTags{curRun,ii}=sprintf('imagename-%s_studyType-%s',...
                            imagename, studyType);
                        
                        clear imagename studyType retType dmScore;
                    catch
                        msgbox('Error specifiying encoding conditions.');
                    end
                    
                case 'TaskB'
                    try
                        %%% pull variables to add to the trial file name
                        % image names
                        
                        % change pic1 to pic for retrieval
                        imageFace = erase(BehavData.pic(currRunIDs(ii)),'.\faces\');
                        imageScene = erase(BehavData.pic2(currRunIDs(ii)),'.\scenes\scenes\');
                        imageFace = erase(imageFace,'"');
                        imageScene = erase(imageScene,'"');
                        imageFace = strrep(imageFace,'''','');
                        imageScene = strrep(imageScene,'''','');
                        
                        % Merge image names
                        imagename = strcat(imageFace{1,1},'|',imageScene{1,1});
                        
                        % Study trial type
                        switch BehavData.familiarnovel(ii)
                            case 0
                                studyType       = 'Trained';
                            case 1
                                studyType       = 'Novel';
                            case 99
                                studyType       = 'Recombined';
                        end
                        
                        % Difference-Memory score (DMscore)
                        switch BehavData.score(ii)
                            case 1 % Miss or False Alarm
                                if BehavData.type(ii) == 0
                                    dmScore = 'Miss';
                                else
                                    dmScore = 'FA';
                                end
                            case 2 % Familiar or Slight False Alarm
                                if BehavData.type(ii) == 0
                                    dmScore = 'Fam-Hit';
                                else
                                    dmScore = 'FA';
                                end
                            case 3 % Hit or Correct Rejection
                                if BehavData.type(ii) == 0
                                    dmScore = 'Rec-Hit';
                                else
                                    dmScore = 'CR';
                                end
                            otherwise
                                dmScore = 'NoBehaviorRecorded';
                        end
                        
                        % informative, unique, BIDS style trial name
                        names{ii} = sprintf...
                            ('imagename-%s_studyType-%s_dmScore%s',...
                            imagename, studyType,  dmScore);
                        
                        ersTags{curRun,ii}=sprintf('imagename-%s_studyType-%s',...
                            imagename, studyType);
                        
                        clear imagename studyType  dmScore;
                    catch
                        msgbox('Error specifiying retrieval conditions.');
                    end
                    
            end
            
        end
        
        %-- Save the Multiple Conditions *.mat file
        % Save the names, onsets, and durations variables to .mat file
        % to be used for later model estimation in SPM. See EstimateModel.m
        
        matfilename = fullfile(curSubj.directory, ...
            ['Run', num2str(curRun, '%03d'), '_multiple_conditions.mat']);
        fprintf('Saving subject %s run %d multiple conditions file...\n',...
            Subjects{i}, curRun)
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
            finalSummary{i+1,1}=Subjects{i};
            finalIndex=[2:4];
        end
        
        finalSummary{i+1,finalIndex(1)}=num2str(length(names));
        finalSummary{i+1,finalIndex(2)}=num2str(length(onsets));
        finalSummary{i+1,finalIndex(3)}=num2str(length(durations));
        
        finalIndex=finalIndex+3;
        
    end
    
    ersTagFilename = fullfile(curSubj.directory,'ersTags.mat');
    save(ersTagFilename,'ersTags');
        
end


%% Write output summary file

file = fopen([study_path filesep 'allTrialsSummary.csv'], 'w');

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
