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

%% Pre-Analysis Setup

% Add CoSMoMVPA to the MATLAB search path
%addpath(genpath('/path/to/CoSMoToolbox'));
addpath(genpath...
    ('/gpfs/group/nad12/default/nad12/toolbox/CoSMoMVPA-master/CoSMoMVPA-master'));

% turn cosmo warnings off
cosmo_warning('off');

%% Set Analysis Parameters & Paths
% Load subject IDs, ROIs, and Condition flags
if exist('flag','var')==0
    
    %Select parameter file is flag does not exist
    uiopen('*.mat')
    
end

% path to save results into
parentDir = fileparts(study_path);
if ~exist('subConds','var')
    switch analysisType
        case 'Searchlight'
            analysis = strcat...
                (analysisName,'_Searchlight_Results_',conds{1,1},'_',conds{1,2});
        otherwise
            analysis = strcat...
                (analysisName,'_Results_',conds{1,1},'_',conds{1,2});
    end
else
    switch analysisType
        case 'Searchlight'
            analysis = strcat...
                (analysisName,'_Searchlight_Results_',conds{1,1},'_',...
                conds{1,2},'_',subConds{1,1});
        otherwise
            analysis = strcat...
                (analysisName,'_Results_',conds{1,1},'_',...
                conds{1,2},'_',subConds{1,1});
    end
end
out_path  = fullfile(parentDir, analysis);

%% Main Body
for iteration=1:length(subjects)*length(rois)
    
    % Loop count for all subject/region combinations
    if iteration==1
        subjCount=1;
        regionCount=1;
    elseif mod(iteration,length(rois))==mod(1,length(rois))
        subjCount=subjCount+1;
        regionCount=1;
    elseif mod(iteration,length(rois))>=0
        regionCount=regionCount+1;
    end
    
    subject = subjects{subjCount};
    ROI     = rois{regionCount};
    
    %% Subject-Specific Directories
    
    % Current subject data paths:
    %  data_path = fullpath to this subject's Single Trial Model directory
    %  spm_path  = fullpath to this subject's SPM.mat file. Note: the
    %                 :beta appended to the end tells cosmo to pull the beta
    %                 information from the SPM.mat file.
    data_path   = fullfile(study_path, subject);
    output_path = fullfile([out_path '_' classificationType], subject);
    spm_fn = [data_path '/SPM.mat'];
    
    % create the output path if it doesn't already exist
    if ~exist(output_path, 'dir')
        mkdir(output_path)
    end
    
    % Path to current region mask
    curROI = fullfile(roi_path, ROI);
    %curROI = fullfile(roi_path, subject, ROI);
    
    switch classificationType
        case 'RSA'
            %%% Loading ROI data into CosmoMVPA - use SPM betas
            % Note that loading data through the SPM.mat file will automatically
            % "chunk" by runs, which is what we want
            fprintf('Loading data from ROI: %s\n',ROI);
            currDataset=cosmo_fmri_dataset([spm_fn ':beta'],'mask',curROI);
            
            %%% Tidy up the dataset
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
    switch classificationType
        case 'RSA'
            try
                if exist('subConds','var')
                    
                    % Comment out for non-ICEE
                    subConds{1,2}=[];
                    subConds={subConds{1,1}};
                    
                    for ii=1:length(conds)
                        
                        %Orig
                        %subCond.(conds{1,ii})=contains(currDataset.sa.labels, conds{1,ii});
                        
                        %ICEE
                        subCond.(conds{1,ii})=contains(currDataset.sa.labels, ['-' conds{1,ii}]);
                        counter=1;
                        
                        for iii=1:length(subConds)
                            subCond.(conds{1,ii})(:,counter+1)=contains...
                                (currDataset.sa.labels, subConds{1,iii});
                            counter=counter+1;
                        end
                        
                        subCond.(conds{1,ii})=double(subCond.(conds{1,ii}));
                        subCond.(conds{1,ii})(:,counter+1)=sum(subCond.(conds{1,ii})(:,1:counter),2);
                        
                    end
                    
                    Cond(1).idx = find(subCond.(conds{1,1})(:,counter+1) == counter);
                    Cond(2).idx = find(subCond.(conds{1,2})(:,counter+1) == counter);
                    
                    CondList = zeros(size(currDataset.samples,1),1);
                    for ii=1:length(conds)
                        
                        CondList(Cond(ii).idx) = ii;
                        
                    end
                    
                else
                    
                    CondList = zeros(size(currDataset.samples,1),1);
                    for ii=1:length(conds)
                        Cond(ii).labels = ~cellfun(@isempty, strfind...
                            (currDataset.sa.labels, conds{ii}));
                        Cond(ii).idx = find(Cond(ii).labels == 1);
                        
                        CondList(Cond(ii).idx) = ii;
                        
                    end
                    
                end
                
                currDataset.sa.targets = CondList;
                
                
                %Codes trials/conditions of no interest as 0 (see SpecifyModel script
                %for trial tag information)
                Zeroidx = find(CondList == 0);
                
                %Removes all trials of no interest from analysis
                if isempty(Zeroidx)==0
                    currDataset.samples(Zeroidx,:)=[];
                    currDataset.sa.targets(Zeroidx)=[];
                    currDataset.sa.beta_index(Zeroidx)=[];
                    currDataset.sa.chunks(Zeroidx)=[];
                    currDataset.sa.fname(Zeroidx)=[];
                    currDataset.sa.labels(Zeroidx)=[];
                end
            catch
                warning('Failure to define trial information. Set to debug mode.');
            end
            
            fprintf('Number of possible targets: %i\n',...
                length(unique(currDataset.sa.targets)));
            
            % Print dataset
            fprintf('Dataset input:\n');
            
            % Print dataset
            fprintf('Number of samples: %i\n',...
                size(currDataset.samples,1));
            fprintf('Number of features (voxels): %i\n',...
                size(currDataset.samples,2));
            fprintf('Number of chunks (runs): %i\n',...
                length(unique(currDataset.sa.chunks)));
    end
    
    %% RSA ROI analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    switch classificationType
        % Single task RSA
        case 'RSA'
            switch analysisType
                case 'Searchlight'
                    
                    %COMING SOON!!!!
                    
                case 'ROI'
                    
                    % Calculate separate matrices for each condition & compare
                    % those
                    
                    switch trialAnalysis
                        case 'Individual'
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
                            
                            % Setup DSM and searchlight arguments
                            % DSM
                            dsmArgs.metric = 'correlation';
                            dsmArgs.center_data = 1;
                            
                            % Create Target DSM from Condtion 1
                            targetDSM = cosmo_dissimilarity_matrix_measure(Conditions(1), dsmArgs);
                            
                            % Set target DSM
                            dsmArgs.target_dsm = targetDSM.samples;
                            
                            % Searchlight ERS for each condition separately
                            rho = cosmo_target_dsm_corr_measure(Conditions(2),dsmArgs);
                            
                            clear targetDSM Conditions;
                            
                        case 'Mean'
                            %args.center_data=true;
                            args.metric='correlation';
                            args.type = 'Pearson';
                            
                            currDataset_mean = cosmo_fx...
                                (currDataset,@(x)mean(x,1),'targets');
                            currDataset_mean.samples(1,:) = zscore(currDataset_mean.samples(1,:));
                            currDataset_mean.samples(2,:) = zscore(currDataset_mean.samples(2,:));
                            currDatasetDSM = cosmo_dissimilarity_matrix_measure...
                                (currDataset_mean,args);
                            
                            rho = (cosmo_squareform(currDatasetDSM.samples) - 1) * -1;
                    end
            end
            
            % ERS analysis
        case 'ERS'
            switch analysisType
                case 'Searchlight'
                    %%% Loading ROI data into CosmoMVPA - use SPM betas
                    fprintf('Loading data from ROI: %s\n',ROI);
                    for i=1:length(tasks)
                        
                        spm_fn = [parentDir '/SingleTrialModel' tasks{i} '/' ...
                            subject '/SPM.mat'];
                        currDataset(i)=cosmo_fmri_dataset([spm_fn ':beta'],'mask',curROI);
                        
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
                    
                    for i=1:length(tasks)
                        currDataset(i).samples(:,logical(removeVoxels))=[];
                        currDataset(i).fa.i(:,logical(removeVoxels))=[];
                        currDataset(i).fa.j(:,logical(removeVoxels))=[];
                        currDataset(i).fa.k(:,logical(removeVoxels))=[];
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
                            
                            if i==1
                                retData(ii)=currDataset(i);
                                retData(ii).samples=retData(ii).samples(Cond(ii).idx,:);
                                retData(ii).sa.beta_index=retData(ii).sa.beta_index(Cond(ii).idx);
                                retData(ii).sa.chunks=retData(ii).sa.chunks(Cond(ii).idx);
                                retData(ii).sa.fname=retData(ii).sa.fname(Cond(ii).idx);
                                retData(ii).sa.labels=retData(ii).sa.labels(Cond(ii).idx);
                                retData(ii).sa.targets=[1:length(Cond(ii).idx)]';
                                
                            else
                                encData(ii)=currDataset(i);
                                encData(ii).samples=encData(ii).samples(Cond(ii).idx,:);
                                encData(ii).sa.beta_index=encData(ii).sa.beta_index(Cond(ii).idx);
                                encData(ii).sa.chunks=encData(ii).sa.chunks(Cond(ii).idx);
                                encData(ii).sa.fname=encData(ii).sa.fname(Cond(ii).idx);
                                encData(ii).sa.labels=encData(ii).sa.labels(Cond(ii).idx);
                                encData(ii).sa.targets=[1:length(Cond(ii).idx)]';
                                
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
                    
                    save([output_path '/searchlightResults_' metric '_' ...
                        num2str(searchlightSize) '.mat'],'results');
                    
                    for i=1:length(conds)
                        % Define output location
                        if ~exist('subConds','var')
                            output{i}=strcat(output_path,'/',subject,...
                                '_ERS_Searchlight_',metric,'_',...
                                num2str(searchlightSize),'_',...
                                conds{i},'.nii');
                        else
                            output=strcat(output_path,'/',subject,'_',...
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
                        
                        spm_fn = [parentDir '/SingleTrialModel' tasks{i} '/' ...
                            subject '/SPM.mat'];
                        currDataset(i)=cosmo_fmri_dataset([spm_fn ':beta'],'mask',curROI);
                        
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
                                    Cond(ii).(tasks{i}).labels = ~cellfun(@isempty, strfind...
                                        (currDataset(i).sa.labels, ['-' conds{ii}]));
                                    Cond(ii).(tasks{i}).idx = find(Cond(ii).(tasks{i}).labels == 1);
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
                                case 'Enc'
                                    encData(ii) = currDataset(i);
                                    encData(ii).samples=encData(ii).samples(Cond(ii).(tasks{i}).idx,:);
                                    encData(ii).sa.beta_index=encData(ii).sa.beta_index(Cond(ii).(tasks{i}).idx);
                                    encData(ii).sa.chunks=encData(ii).sa.chunks(Cond(ii).(tasks{i}).idx);
                                    encData(ii).sa.fname=encData(ii).sa.fname(Cond(ii).(tasks{i}).idx);
                                    encData(ii).sa.labels=encData(ii).sa.labels(Cond(ii).(tasks{i}).idx);
                                    encData(ii).sa.targets=[1:length(Cond(ii).(tasks{i}).idx)]';
                                case 'Ret'
                                    retData(ii) = currDataset(i);
                                    retData(ii).samples = retData(ii).samples(Cond(ii).(tasks{i}).idx,:);
                                    retData(ii).sa.beta_index = retData(ii).sa.beta_index(Cond(ii).(tasks{i}).idx);
                                    retData(ii).sa.chunks = retData(ii).sa.chunks(Cond(ii).(tasks{i}).idx);
                                    retData(ii).sa.fname = retData(ii).sa.fname(Cond(ii).(tasks{i}).idx);
                                    retData(ii).sa.labels = retData(ii).sa.labels(Cond(ii).(tasks{i}).idx);
                                    retData(ii).sa.targets = [1:length(Cond(ii).(tasks{i}).idx)]';
                            end
                        end
                    end
                    
                    % Setup DSM and searchlight arguments
                    %DSM
                    dsmArgs.metric = 'correlation';
                    dsmArgs.center_data = 1;
                    
                    for i=1:length(Cond)
                        
                        % Create Target DSM from Encoding Run
                        targetDSM(i) = cosmo_dissimilarity_matrix_measure(encData(i), dsmArgs);
                        
                        % Set target DSM
                        dsmArgs.target_dsm = targetDSM(i).samples;
                        
                        % Searchlight ERS for each condition separately
                        rho.(conds{i})  = cosmo_target_dsm_corr_measure(retData(i),dsmArgs);
                        
                    end
                    
                    clear remove removeVoxels encData retData;
                    
            end
    end
    
    %% Save text output of SVM Classification
    if strcmpi(analysisType,'Searchlight')==0
        % Create a tidyverse formatted table for final statistical analysis
        TrialTypeCombo = {strcat(conds{1,1},'_v_',conds{1,2})};
        
        % create subjectid and roiid columns
        regionName=erase(ROI,{'reslice_','_bilat.nii'});
        subjectid   = repmat({subject}, length(TrialTypeCombo), 1);
        roiid       = repmat({regionName}, length(TrialTypeCombo), 1);
        
        switch classificationType
            case 'RSA'
                % create the stats table
                %stats_table = table...
                %    (subjectid, roiid, TrialTypeCombo, rho(1,2));
                stats_table = table...
                    (subjectid, roiid, TrialTypeCombo, rho.samples);
            case 'ERS'
                % create the stats table
                stats_table = table...
                    (subjectid, roiid, TrialTypeCombo, ...
                    rho.(conds{1}).samples, rho.(conds{2}).samples);
                %rho.(conds{1})(1,2), rho.(conds{2})(1,2));
        end
        
        % write the stats table
        filename = sprintf('sub-%s_roiid-%s_statistics-table.csv', subject, ROI);
        writetable(stats_table, fullfile(output_path, filename));
        
        %% Create aggregate table for easy viewing
        % Create headers
        if iteration==1
            
            switch classificationType
                case 'RSA'
                    finalTable=cell(length(subjects)+1,length(rois)+2);
                    finalTable{1,1}='subjectid';
                    finalTable{1,2}='Trial Type';
                    tempcount=3;
                    
                    for header=1:length(rois)
                        finalTable{1,tempcount}=strcat(...
                            rois{1,header}(1:end-4),'_Similarity');
                        tempcount=tempcount+1;
                    end
                case 'ERS'
                    finalTable=cell(length(subjects)+1,length(rois)*2+2);
                    finalTable{1,1}='subjectid';
                    finalTable{1,2}='Trial Type';
                    tempcount=3;
                    
                    for header=1:length(rois)
                        finalTable{1,tempcount}=strcat(...
                            rois{1,header}(1:end-4),'_',conds{1,1},'_Dissimilarity');
                        finalTable{1,tempcount+1}=strcat(...
                            rois{1,header}(1:end-4),'_',conds{1,2},'_Dissimilarity');
                        tempcount=tempcount+2;
                    end
            end
            
            finalTable{1,tempcount}=['num' conds{1,1}];
            finalTable{1,tempcount+1}=['num' conds{1,2}];
            
            row=2;
            header=3;
            clear tempcount;
        end
        
        % Counter for resetting to next row. Uses remainder from divison of
        % region counter over total (e.g. 1/14) to check data should be
        % read into next subject line/row.
        iterCheck=mod(iteration,length(rois));
        
        % Add subject, trial type, and accuracy to table
        finalTable{row,1}=subject;
        finalTable{row,2}=TrialTypeCombo{1,1};
        switch classificationType
            case 'RSA'
                %finalTable{row,header}=currDatasetDSM.samples;
                finalTable{row,header}=rho.samples;
                header=header+1;
            case 'ERS'
                %finalTable{row,header}=rho.(conds{1})(1,2);
                %finalTable{row,header+1}=rho.(conds{2})(1,2);
                
                finalTable{row,header}=rho.(conds{1}).samples;
                finalTable{row,header+1}=rho.(conds{2}).samples;
                
                header=header+2;
        end
        
        % Drops to next row if remainder is 0 (e.g. all regions have been
        % entered for a given subject)
        if iterCheck == 0
            finalTable{row,header}=num2str(length(Cond(1).(tasks{i}).idx));
            finalTable{row,header+1}=num2str(length(Cond(2).(tasks{i}).idx));
            row=row+1;
            header=3;
        end
        
    end
end

% Save mat file with statistics separately
save([out_path '_' classificationType filesep 'finalTable.mat'],'finalTable');

%% Save summary files of RSA.
if strcmpi(analysisType,'Searchlight')==0
    
    % Write output summary file
    file = fopen([out_path '_' classificationType filesep 'allAccuraciesSummary.csv'], 'w');
    
    for a=1:size(finalTable,1)
        for b=1:size(finalTable,2)
            var = eval('finalTable{a,b}');
            try
                fprintf(file, '%s', var);
            end
            fprintf(file, ',');
        end
        fprintf(file, '\n');
    end
    
    fclose(file);
    clear;
end
