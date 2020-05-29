%% Run Representational Similarity Analysis (RSA)
%   Editor:    Daniel Elbich
%   Updated:   4/29/19
%
% Representational similarity analysis (RSA) for a single subject. Flagged
% for either single ROI SVM classification or searchlight analysis.
%
% Load single-trial beta images from each subject, apply ROI mask, calculate
% correlations between trial patterns, take the mean across trial types
%
% Current Developer: Daniel Elbich, delbich10@gmail.com
% 2/22/19
%
% Updates:
%
% 3/28/19 - Loads in parameter file created by createParams.m subscript. Mat
% file should contain paths to roi and data folders, lists of rois,
% subjects, & conditions, and analysis name. See createParams.m for full
% list.

%% Set Analysis Parameters & Paths
% Load all relevent project information
if exist('commandFlag','var') == 0
    
    %Select parameter file is flag does not exist
    [file,path]=uigetfile('*.mat','Select params file');
    filename=fullfile(path,file);
    load(filename);
    
end

% turn cosmo warnings off
cosmo_warning('off');

% Filepath for results folder
parentDir = directory.Model;

% Base output directory name
analysis = [directory.Analysis filesep 'models' filesep file(1:end-4)];

%Debug
%subjects(2)=[];

%% Main Body
for iteration=1:length(subjects)
    
    %% Subject-Specific Directories
    
    % Current subject data paths:
    %  dataPath = fullpath to this subject's Single Trial Model directory
    %  spmFile  = fullpath to this subject's SPM.mat file. Note: the
    %                 :beta appended to the end tells cosmo to pull the beta
    %                 information from the SPM.mat file.
    dataPath   = fullfile(directory.Model, subjects{iteration});
    outputPath = fullfile(analysis, subjects{iteration});
    spmFile = [dataPath '/SPM_gz.mat'];
    
    % create the output path if it doesn't already exist
    if ~exist(outputPath, 'dir')
        mkdir(outputPath)
    end
    
    %% Load Mask Data
    masks = dir([directory.Analysis filesep 'masks' filesep file(1:end-4)...
        filesep subjects{iteration} filesep '*.nii.gz']);
    
    %Debug
    masks=masks(1:3);
    
    for curMask = 1:length(masks)
        
        switch classType
            case 'RSA'
                % Path to current region mask
                curROI = fullfile(masks(curMask).folder, masks(curMask).name);
                
                % Current region name
                regionName=erase(masks(curMask).name,'.nii.gz');
                
                %%% Loading ROI data into CosmoMVPA - use SPM betas
                % Note that loading data through the SPM.mat file will automatically
                % "chunk" by runs, which is what we want
                fprintf('Loading data from ROI: %s\n',masks(curMask).name);
                currDataset=cosmo_fmri_dataset([spmFile ':beta'],'mask',curROI);
                
                % Clear errant Not a Number (NaN) values
                % Remove constant features
                currDataset=cosmo_remove_useless_data(currDataset);
                
                switch regressRT.flag
                    case 'Yes'
                        files = dir([study_path filesep subject filesep 'Run*']);
                        for i=1:length(files)
                            curMat(i) = load([files(i).folder filesep files(i).name]);
                            if i==length(files)
                                rtCell = [curMat.RT];
                                
                                % Convert from cell to double for regression
                                for ii=1:length(rtCell)
                                    
                                    % Flag outlier RT greater than 4 seconds
                                    if double(rtCell{ii}) >= regressRT.trialSec
                                        rtDouble(ii,1) = regressRT.trialSec;
                                    else
                                        rtDouble(ii,1) = double(rtCell{ii});
                                    end
                                end
                                
                                % Replace with trial duration (set in params)
                                rtDouble(isnan(rtDouble))=regressRT.trialSec;
                            end
                        end
                        
                        for i=1:length(currDataset.samples)
                            model = LinearModel.fit(rtDouble,currDataset.samples(:,i));
                            if i==1
                                allResiduals = model.Residuals.Raw;
                            else
                                allResiduals = [allResiduals model.Residuals.Raw];
                            end
                        end
                        
                        zscoreResid = zscore(allResiduals);
                        currDataset.samples = zscoreResid;
                        
                        clear files curMat rtCell rtDouble model allResiduals zscoreResid;
                end
            case 'ERS'
                disp('Skipping to ERS...');
        end
        
        %% Define trial information
        switch classType
            case 'RSA'
                try
                    % Identify trials of interest for each condition
                    if exist('subConds','var')
                        
                        for ii=1:length(taskInfo.Conditions)
                            
                            subCond.(taskInfo.Conditions{1,ii})=contains(currDataset.sa.labels, taskInfo.Conditions{1,ii});
                            counter=1;
                            
                            for iii=1:length(subConds)
                                subCond.(taskInfo.Conditions{1,ii})(:,counter+1)=contains...
                                    (currDataset.sa.labels, subConds{1,iii});
                                counter=counter+1;
                            end
                            
                            subCond.(taskInfo.Conditions{1,ii})=double(subCond.(taskInfo.Conditions{1,ii}));
                            subCond.(taskInfo.Conditions{1,ii})(:,counter+1)=sum(subCond.(taskInfo.Conditions{1,ii})(:,1:counter),2);
                            
                        end
                        
                        Cond(1).idx = find(subCond.(taskInfo.Conditions{1,1})(:,counter+1) == counter);
                        Cond(2).idx = find(subCond.(taskInfo.Conditions{1,2})(:,counter+1) == counter);
                        
                    else
                        
                        CondList = zeros(size(currDataset.samples,1),1);
                        for ii=1:length(taskInfo.Conditions)
                            Cond(ii).labels = ~cellfun(@isempty, strfind...
                                (currDataset.sa.labels, taskInfo.Conditions{ii}));
                            Cond(ii).idx = find(Cond(ii).labels == 1);
                            
                            CondList(Cond(ii).idx) = ii;
                            
                        end
                        
                    end
                    
                    currDataset.sa.targets = CondList;
                    
                    % Codes trials/conditions of no interest as 0 (see SpecifyModel script
                    % for trial tag information)
                    Zeroidx = find(CondList == 0);
                    
                    % Removes all trials of no interest from analysis
                    if isempty(Zeroidx)==0
                        currDataset.samples(Zeroidx,:)=[];
                        currDataset.sa.targets(Zeroidx)=[];
                        currDataset.sa.beta_index(Zeroidx)=[];
                        currDataset.sa.chunks(Zeroidx)=[];
                        currDataset.sa.fname(Zeroidx)=[];
                        currDataset.sa.labels(Zeroidx)=[];
                    end
                    
                    fprintf('Number of possible targets: %i\n',length...
                        (unique(currDataset.sa.targets)));
                    
                    % Print dataset
                    fprintf('Dataset input:\n');
                    
                    % Print dataset
                    fprintf('Number of samples: %i\n',...
                        size(currDataset.samples,1));
                    fprintf('Number of features (voxels): %i\n',...
                        size(currDataset.samples,2));
                    fprintf('Number of chunks (runs): %i\n',...
                        length(unique(currDataset.sa.chunks)));
                    
                    % RSA ROI analysis
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    switch analysisType
                        case 'Searchlight'
                            
                            %%% IN PROGRESS %%%
                            
                        case 'ROI'
                            
                            % Calculate separate matrices for each condition & compare
                            % those
                            
                            % Split dataset into separate condition
                            % variables
                            for i=1:length(Cond)
                                index = find(currDataset.sa.targets == i);
                                Conditions(i) = currDataset;
                                Conditions(i).samples = Conditions(i).samples(index,:);
                                Conditions(i).sa.beta_index = Conditions(i).sa.beta_index(index);
                                Conditions(i).sa.chunks = Conditions(i).sa.chunks(index);
                                Conditions(i).sa.fname = Conditions(i).sa.fname(index);
                                Conditions(i).sa.labels = Conditions(i).sa.labels(index);
                                
                                % Re-index trials within condition into
                                % separate targets
                                Conditions(i).sa.targets = [1:length(index)]';
                            end
                            
                            % Re-index trials within condition into
                            % separate targets
                            Conditions(i).sa.targets = [1:length(index)]';
                            
                            % Setup DSM and searchlight arguments
                            % DSM
                            
                            dsmArgs(i).metric = 'correlation';
                            dsmArgs(i).center_data = 1;
                            dsmArgs(i).trialType = taskInfo.Conditions(1,i);
                            
                            % Create Target DSM from Condtion 1
                            targetDSM(i) = cosmo_dissimilarity_matrix_measure(Conditions(i), dsmArgs(i));
                            
                            %Convert to similarity for all conds
                            %individualls
                            for ii=1:length(targetDSM(i).samples)
                                if targetDSM(i).samples(ii)<1 && targetDSM(i).samples(ii)>0
                                    newSim(ii,i) = 1 - targetDSM(i).samples(ii);
                                else
                                    newSim(ii,i) = -targetDSM(i).samples(ii);
                                    newSim(ii,i) = newSim(ii,i) + 1;
                                end
                            end
                            
                            % Set target DSM
                            dsmArgs(i).target_dsm = targetDSM(i).samples;
                    end
                    
                    %New
                    newRow=1;
                    for i=1:length(taskInfo.Conditions)
                        meanCorrWithinCondition{newRow,i} = mean(newSim(:,i));
                    end
                    newRow = newRow + 1;
                    
                    %end new
                    
                    % Get raw correlation matrix back
                    % targetDSMCorrMatrix = cosmo_squareform(targetDSM.samples)
                    % newCorrWithinMatrix = cosmo_squareform(newSim)
                    
                    % Obtain number of combinations with 2
                    % conditions selected
                    combinations=factorial(length(Cond))/...
                        ((factorial(2)*(factorial(length(Cond)-2))));
                    
                    % Searchlight ERS for each condition separately
                    iter=1;
                    for i=1:combinations
                        if i<length(Cond)
                            condCount=1;
                        elseif i<combinations
                            condCount=2;
                            if i==length(Cond)
                                iter=condCount;
                            end
                        else
                            condCount=3;
                            iter=condCount;
                        end
                        
                        rho(1,i) = cosmo_target_dsm_corr_measure...
                            (Conditions(condCount),dsmArgs(iter+1));
                        
                        TrialTypeCombo(1,i) = {strcat(taskInfo.Conditions{condCount},'_v_',taskInfo.Conditions{iter+1})};
                        
                        iter=iter+1;
                    end
                    
                    
                    clear targetDSM Conditions;
                end
                
                %% ERS analysis
            case 'ERS'
                switch analysisType
                    case 'Searchlight'
                        %%% Loading ROI data into CosmoMVPA - use SPM betas
                        fprintf('Loading data from ROI: %s\n',ROI);
                        for i=1:length(tasks)
                            
                            spmFile = [parentDir '/SingleTrialModel' tasks{i} '/' ...
                                subject '/SPM_gz.mat'];
                            spmFile = [dataPath '/SPM_gz.mat'];
                            currDataset(i)=cosmo_fmri_dataset([spmFile ':beta'],'mask',curROI);
                            
                        end
                        
                        % Remove NaN voxels in both datasets
                        for i=1:size(currDataset(1).samples,2)
                            for j=1:length(tasks)
                                remove(j) = isnan(currDataset(j).samples(1,i));
                            end
                            
                            if sum(remove)>=1
                                removeVoxels(i) = 1;
                            else
                                removeVoxels(i) = 0;
                            end
                        end
                        
                        for i=1:length(tasks)
                            currDataset(i).samples(:,logical(removeVoxels))=[];
                            currDataset(i).fa.i(:,logical(removeVoxels))=[];
                            currDataset(i).fa.j(:,logical(removeVoxels))=[];
                            currDataset(i).fa.k(:,logical(removeVoxels))=[];
                            
                            if strcmpi(task{i},'Retrieval')==1
                                switch regressRT.flag
                                    case 'Yes'
                                        files = dir([study_path filesep subject filesep 'Run*']);
                                        for j=1:length(files)
                                            curMat(j) = load([files(j).folder filesep files(j).name]);
                                            if j==length(files)
                                                rtCell = [curMat.RT];
                                                
                                                % Convert from cell to double for regression
                                                for k=1:length(rtCell)
                                                    
                                                    % Flag outlier RT greater than 4 seconds
                                                    if double(rtCell{k}) >= regressRT.trialSec
                                                        rtDouble(k,1) = regressRT.trialSec;
                                                    else
                                                        rtDouble(k,1) = double(rtCell{k});
                                                    end
                                                end
                                                
                                                % Replace with trial duration (set in params)
                                                rtDouble(isnan(rtDouble))=regressRT.trialSec;
                                            end
                                        end
                                        
                                        % Z-Score RT values
                                        rtDouble = zscore(rtDouble);
                                        
                                        for jj=1:length(currDataset.samples)
                                            model = LinearModel.fit(rtDouble,currDataset.samples(:,j));
                                            if j==1
                                                allResiduals = model.Residuals.Raw;
                                            else
                                                allResiduals = [allResiduals model.Residuals.Raw];
                                            end
                                        end
                                        
                                        zscoreResid = zscore(allResiduals);
                                        currDataset.samples = zscoreResid;
                                        
                                        clear files curMat rtCell rtDouble model allResiduals zscoreResid;
                                end
                                
                            end
                            
                            try
                                if exist('subConds','var')
                                    
                                    subConds{1,2}=[];
                                    subConds={subConds{1,1}};
                                    
                                    for ii=1:length(conds)
                                        
                                        subCond(ii).idx=contains(currDataset(i).sa.labels, ['-' conds{ii}]);
                                        counter=1;
                                        
                                        for iii=1:length(subConds)
                                            subCond(ii).idx(:,counter+1)=contains(currDataset(i).sa.labels, subConds{1,iii});
                                            counter=counter+1;
                                        end
                                        
                                        subCond(ii).idx=double(subCond(ii).idx);
                                        subCond(ii).idx(:,counter+1)=sum(subCond(ii).idx(:,1:counter),2);
                                        
                                    end
                                    
                                    Cond(1).idx = find(subCond(1).idx(:,counter+1) == counter);
                                    Cond(2).idx = find(subCond(2).idx(:,counter+1) == counter);
                                    
                                else
                                    
                                    for ii=1:length(conds)
                                        Cond(ii).labels = ~cellfun(@isempty, strfind...
                                            (currDataset(i).sa.labels, ['-' conds{ii}]));
                                        Cond(ii).idx = find(Cond(ii).labels == 1);
                                    end
                                    
                                end
                                
                            catch
                                warning('Failure to define trial information. Set to debug mode.');
                            end
                        end
                        
                        % Parse data by task and condition
                        counter=1;
                        for i=1:length(currDataset)
                            for ii=1:length(Cond)
                                if strcmpi(task{i},'Retrieval')==1
                                    retData(ii)=currDataset(i);
                                    retData(ii).samples=retData(ii).samples(Cond(ii).idx,:);
                                    retData(ii).sa.beta_index=retData(ii).sa.beta_index(Cond(ii).idx);
                                    retData(ii).sa.chunks=retData(ii).sa.chunks(Cond(ii).idx);
                                    retData(ii).sa.fname=retData(ii).sa.fname(Cond(ii).idx);
                                    retData(ii).sa.labels=retData(ii).sa.labels(Cond(ii).idx);
                                    retData(ii).sa.targets=retData(ii).sa.targets(Cond(ii).idx);
                                    %retData(ii).sa.targets=[1:length(Cond(ii).idx)]';
                                else
                                    encData(ii)=currDataset(i);
                                    encData(ii).samples=encData(ii).samples(Cond(ii).idx,:);
                                    encData(ii).sa.beta_index=encData(ii).sa.beta_index(Cond(ii).idx);
                                    encData(ii).sa.chunks=encData(ii).sa.chunks(Cond(ii).idx);
                                    encData(ii).sa.fname=encData(ii).sa.fname(Cond(ii).idx);
                                    encData(ii).sa.labels=encData(ii).sa.labels(Cond(ii).idx);
                                    encData(ii).sa.targets=encData(ii).sa.targets(Cond(ii).idx);
                                    %encData(ii).sa.targets=[1:length(Cond(ii).idx)]';
                                end
                            end
                        end
                        
                        % Setup DSM and searchlight arguments
                        %DSM
                        dsmArgs.metric = 'correlation';
                        dsmArgs.center_data = 1;
                        
                        %Searchlight
                        searchlightSize = 20;
                        metric = 'count';
                        measure = @cosmo_target_dsm_corr_measure;
                        measure_args = struct();
                        measure_args.type = 'Spearman';
                        measure_args.center_data = 1;
                        
                        
                        for i=1:length(Cond)
                            
                            % Create Target DSM from Encoding Run
                            targetDSM(i) = cosmo_dissimilarity_matrix_measure(encData(i), dsmArgs);
                            
                            if i==1
                                % Create Searchlight
                                nbrhood=cosmo_spherical_neighborhood(retData(i),metric,searchlightSize);
                            end
                            
                            % Set target DSM
                            measure_args.target_dsm = targetDSM(i).samples;
                            
                            % Searchlight ERS for each condition separately
                            results(i) = cosmo_searchlight(retData(i),nbrhood,measure,measure_args);
                            
                        end
                        
                        save([outputPath '/searchlightResults_' metric '_' ...
                            num2str(searchlightSize) '.mat'],'results');
                        
                        for i=1:length(conds)
                            % Define output location
                            if ~exist('subConds','var')
                                output{i}=strcat(outputPath,'/',subject,...
                                    '_ERS_Searchlight_',metric,'_',...
                                    num2str(searchlightSize),'_',...
                                    conds{i},'.nii');
                            else
                                output=strcat(outputPath,'/',subject,'_',...
                                    classifier.name,'_Searchlight_',metric,'_',...
                                    num2str(searchlightSize),'_',...
                                    conds{1,1},'_vs_',conds{1,2},'_',subConds{1,1},'.nii');
                            end
                            
                            % Store results to disk
                            cosmo_map2fmri(results(i), output{i});
                            
                        end
                        
                    case 'ROI'
                        %%% Loading ROI data into CosmoMVPA - use SPM betas
                        fprintf('Loading data from ROI: %s\n',ROI);
                        for i=1:length(tasks)
                            
                            spmFile = [parentDir '/SingleTrialModel' tasks{i} '/' ...
                                subject '/SPM.mat'];
                            currDataset(i)=cosmo_fmri_dataset([spmFile ':beta'],'mask',curROI);
                            
                        end
                        
                        % Remove errant voxels in both datasets
                        for i=1:size(currDataset(1).samples,2)
                            for j=1:length(tasks)
                                remove(j) = isnan(currDataset(j).samples(1,i));
                            end
                            
                            if sum(remove)>=1
                                removeVoxels(i) = 1;
                            else
                                removeVoxels(i) = 0;
                            end
                        end
                        
                        %%% Indentify trials of interest for each condition
                        for i=1:length(tasks)
                            currDataset(i).samples(:,logical(removeVoxels))=[];
                            currDataset(i).fa.i(:,logical(removeVoxels))=[];
                            currDataset(i).fa.j(:,logical(removeVoxels))=[];
                            currDataset(i).fa.k(:,logical(removeVoxels))=[];
                            
                            % Mean center data
                            currDataset(i).samples = bsxfun...
                                (@minus,currDataset(i).samples,mean(currDataset(i).samples,1));
                            
                            if strcmpi(tasks{i},'Retrieval')==1
                                switch regressRT.flag
                                    case 'Yes'
                                        files = dir([study_path filesep subject filesep 'Run*']);
                                        for j=1:length(files)
                                            curMat(j) = load([files(j).folder filesep files(j).name]);
                                            if j==length(files)
                                                rtCell = [curMat.RT];
                                                
                                                % Convert from cell to double for regression
                                                for k=1:length(rtCell)
                                                    
                                                    % Flag outlier RT greater than 4 seconds
                                                    if double(rtCell{k}) >= regressRT.trialSec
                                                        rtDouble(k,1) = regressRT.trialSec;
                                                    else
                                                        rtDouble(k,1) = double(rtCell{k});
                                                    end
                                                end
                                                
                                                % Replace with trial duration (set in params)
                                                rtDouble(isnan(rtDouble))=regressRT.trialSec;
                                            end
                                        end
                                        
                                        % Z-Score RT values
                                        rtDouble = zscore(rtDouble);
                                        
                                        for jj=1:length(currDataset.samples)
                                            model = LinearModel.fit(rtDouble,currDataset.samples(:,j));
                                            if j==1
                                                allResiduals = model.Residuals.Raw;
                                            else
                                                allResiduals = [allResiduals model.Residuals.Raw];
                                            end
                                        end
                                        
                                        zscoreResid = zscore(allResiduals);
                                        currDataset.samples = zscoreResid;
                                        
                                        clear files curMat rtCell rtDouble model allResiduals zscoreResid;
                                end
                                
                            end
                            
                            try
                                if exist('subConds','var')
                                    
                                    subConds{1,2}=[];
                                    subConds={subConds{1,1}};
                                    
                                    for ii=1:length(conds)
                                        
                                        subCond(ii).idx=contains(currDataset(i).sa.labels, ['-' conds{ii}]);
                                        counter=1;
                                        
                                        for iii=1:length(subConds)
                                            subCond(ii).idx(:,counter+1)=contains(currDataset(i).sa.labels, subConds{1,iii});
                                            counter=counter+1;
                                        end
                                        
                                        subCond(ii).idx=double(subCond(ii).idx);
                                        subCond(ii).idx(:,counter+1)=sum(subCond(ii).idx(:,1:counter),2);
                                        
                                    end
                                    
                                    Cond(1).idx = find(subCond(1).idx(:,counter+1) == counter);
                                    Cond(2).idx = find(subCond(2).idx(:,counter+1) == counter);
                                    
                                else
                                    
                                    currDataset(i).sa.targets = zeros(size(currDataset(i).samples,1),1);
                                    for ii=1:length(conds)
                                        Cond(ii).(tasks{i}).labels = ~cellfun(@isempty, strfind...
                                            (currDataset(i).sa.labels, conds{ii}));
                                        Cond(ii).(tasks{i}).idx = find(Cond(ii).(tasks{i}).labels == 1);
                                        currDataset(i).sa.targets(Cond(ii).(tasks{i}).idx) = ii;
                                    end
                                    
                                end
                                
                            catch
                                warning('Failure to define trial information. Set to debug mode.');
                            end
                        end
                        
                        % Parse data by task and condition
                        counter=1;
                        for i=1:length(currDataset)
                            for ii=1:length(Cond)
                                switch tasks{i}
                                    case 'Encoding'
                                        encData(ii) = currDataset(i);
                                        encData(ii).samples=encData(ii).samples(Cond(ii).(tasks{i}).idx,:);
                                        encData(ii).sa.beta_index=encData(ii).sa.beta_index(Cond(ii).(tasks{i}).idx);
                                        encData(ii).sa.chunks=encData(ii).sa.chunks(Cond(ii).(tasks{i}).idx);
                                        encData(ii).sa.fname=encData(ii).sa.fname(Cond(ii).(tasks{i}).idx);
                                        encData(ii).sa.labels=encData(ii).sa.labels(Cond(ii).(tasks{i}).idx);
                                        encData(ii).sa.targets=encData(ii).sa.targets(Cond(ii).(tasks{i}).idx);
                                    case 'Retrieval'
                                        retData(ii) = currDataset(i);
                                        retData(ii).samples = retData(ii).samples(Cond(ii).(tasks{i}).idx,:);
                                        retData(ii).sa.beta_index = retData(ii).sa.beta_index(Cond(ii).(tasks{i}).idx);
                                        retData(ii).sa.chunks = retData(ii).sa.chunks(Cond(ii).(tasks{i}).idx);
                                        retData(ii).sa.fname = retData(ii).sa.fname(Cond(ii).(tasks{i}).idx);
                                        retData(ii).sa.labels = retData(ii).sa.labels(Cond(ii).(tasks{i}).idx);
                                        retData(ii).sa.targets = retData(ii).sa.targets(Cond(ii).(tasks{i}).idx);
                                end
                            end
                        end
                        
                        for i=1:length(Cond)
                            
                            for j=1:length(Cond(i).Encoding.idx)
                                for k=1:length(Cond(i).Retrieval.idx)
                                    corrVal(j,k) = corr(encData(i).samples(j,:)',retData(i).samples(k,:)','Type','Pearson');
                                end
                            end
                            
                            rho.(conds{i}) = mean(mean(corrVal));
                            clear corrVal;
                        end
                        
                        clear currDataset removeVoxels encData retData
                end
        end
        
        %% Save text output of SVM Classification
        if strcmpi(analysisType,'Searchlight')==0
            % Create a tidyverse formatted table for final statistical analysis
            
            switch classType
                case 'RSA'
                    % create subjectid and roiid columns
                    subjectid   = repmat(subjects(iteration), length(TrialTypeCombo(1)), 1);
                    roiid       = repmat({regionName}, length(TrialTypeCombo(1)), 1);
                    
                    % create the stats table
                    %stats_table = table...
                    %    (subjectid, roiid, TrialTypeCombo, rho(1,2));
                    
                    for i=1:combinations
                        
                        conditions   = repmat({TrialTypeCombo(i)}, length(TrialTypeCombo(1)), 1);
                        corrVal       = repmat({rho(i).samples}, length(TrialTypeCombo(1)), 1);
                        
                        if i==1
                            statsTable = table(subjectid, roiid, conditions, corrVal);
                            statsTable.Properties.VariableNames{3}=['conditions' num2str(i)];
                            statsTable.Properties.VariableNames{4}=['corrVal' num2str(i)];
                        else
                            tempTable = table(conditions, corrVal);
                            tempTable.Properties.VariableNames{1}=['conditions' num2str(i)];
                            tempTable.Properties.VariableNames{2}=['corrVal' num2str(i)];
                            statsTable = [statsTable tempTable];
                        end
                        
                    end
                    
                    %Add in within-condition correlation
                    statsTable = [statsTable];
                    
                    % write the stats table
                    filename = sprintf('sub-%s_roiid-%s_statistics-table.csv', subjectid{:}, roiid{:});
                    writetable(statsTable, fullfile(outputPath, filename));
                    
                case 'ERS'
                    % Create subjectid and roiid columns
                    TrialTypeCombo = {strcat(tasks{1},'_v_',tasks{2})};
                    regionName  = erase(ROI,{'reslice_','_bilat.nii'});
                    
                    % Create final stats table for region
                    statsTable = cell(2,3+length(conds));
                    for i=1:length(conds)
                        if i==1
                            statsTable(1,1:3) = {'subjectid','roiid','TrialTypeCombo'};
                            statsTable(2,1:3) = {subject, regionName, char(TrialTypeCombo)};
                        end
                        statsTable{1,i+3} = conds{i};
                        statsTable{2,i+3} = rho.(conds{i});
                    end
                    
                    % Write output summary file
                    file = fopen([sprintf([outputPath '/'...
                        'sub-%s_roiid-%s_statistics-table.csv'], subject, regionName)], 'w');
                    
                    for a=1:size(statsTable,1)
                        for b=1:size(statsTable,2)
                            var = eval('statsTable{a,b}');
                            try
                                fprintf(file, '%s', var);
                            end
                            fprintf(file, ',');
                        end
                        fprintf(file, '\n');
                    end
                    fclose(file);
            end
            
            %% Create aggregate table for easy viewing
            % Create headers
            if iteration==1 && curMask==1
                
                switch classType
                    case 'RSA'
                        summary=cell(length(subjects)+1,length(masks)+2);
                        summary{1,1}='subjectid';
                        summary{1,2}='Trial Type';
                        tmpCnt=3;
                        
                        for header=1:length(masks)
                            summary{1,tmpCnt}=[masks(header).name(1:end-7)...
                                '_Similarity'];
                            tmpCnt=tmpCnt+1;
                        end
                    case 'ERS'
                        summary=cell(length(subjects)+1,length(masks)*2+2);
                        summary{1,1}='subjectid';
                        summary{1,2}='Trial Type';
                        tmpCnt=3;
                        
                        for header=1:length(masks)
                            summary{1,tmpCnt}=strcat(masks(header).name(1:end-7),...
                                '_',taskInfo.Conditions{1},'_Similarity');
                            summary{1,tmpCnt+1}=strcat(masks(header).name(1:end-7),...
                                '_',taskInfo.Conditions{2},'_Similarity');
                            tmpCnt=tmpCnt+2;
                        end
                end
                
                summary{1,tmpCnt}=['num_' taskInfo.Conditions{1}];
                summary{1,tmpCnt+1}=['num_' taskInfo.Conditions{2}];
                
                row=2;
                header=3;
                clear tmpCnt;
            end
            
            % Counter for resetting to next row. Uses remainder from divison of
            % region counter over total (e.g. 1/14) to check data should be
            % read into next subject line/row.
            iterCheck=mod(curMask,length(masks));
            
            % Add subject, trial type, and accuracy to table
            summary{row,1}=subjects{iteration};
            summary{row,2}=TrialTypeCombo{1,1};
            switch classType
                case 'RSA'
                    summary{row,header}=rho.samples;
                    header=header+1;
                case 'ERS'
                    summary{row,header}=rho.(taskInfo.Conditions{1});
                    summary{row,header+1}=rho.(taskInfo.Conditions{2});
                    header=header+2;
            end
            
            % Drops to next row if remainder is 0 (e.g. all regions have been
            % entered for a given subject)
            if iterCheck == 0
                switch classType
                    case 'RSA'
                        summary{row,header}=num2str(length(Cond(1).idx));
                        summary{row,header+1}=num2str(length(Cond(2).idx));
                    case 'ERS'
                        summary{row,header}=num2str(length(Cond(1).(tasks{i}).idx));
                        summary{row,header+1}=num2str(length(Cond(2).(tasks{i}).idx));
                end
                row=row+1;
                header=3;
            end
            
        end
    end
end

% Save mat file with statistics separately
save([fileparts(outputPath) filesep 'summary.mat'],'summary');

%% Save summary files of RSA.
if strcmpi(analysisType,'Searchlight')==0
    
    % Write output summary file
    file = fopen([fileparts(outputPath) filesep 'allSimilaritiesSummary.csv'], 'w');
    
    for a=1:size(summary,1)
        for b=1:size(summary,2)
            var = eval('summary{a,b}');
            try
                fprintf(file, '%s', var);
            end
            fprintf(file, ',');
        end
        fprintf(file, '\n');
    end
    
    fclose(file);
    clc;
    clear;
end
