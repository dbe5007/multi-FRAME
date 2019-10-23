%% Run MVPA Classification
%   Editor:    Daniel Elbich
%   Updated:   2/27/19
%
%   Multivoxel pattern analysis (MVPA) for a single subject. Flagged for
%   either single ROI SVM classification or searchlight analysis.
%
%   Load single-trial beta images from each subject, apply ROI mask,
%   calculate correlations between trial patterns, take the mean across
%   trial types.
%
%   Current Developer: Daniel Elbich, delbich10@gmail.com
%   2/22/19
%
%   Updates:
%
%   5/6/19 - Includes flag to collapse across run to include n-1 trials as
%   training set and tests on single trial. Goal is to improve
%   classification by increasing number of training samples.
%   See EstimateModel.m for flag declaration.
%
%   Final predictions by trial by region are saved as mat file for each
%   subject separately. Order will match region testing order.
%
%   3/28/19 - Loads in parameter file created by createParams.m subscript.
%   Mat file should contain paths to roi and data folders, lists of rois,
%   subjects, & conditions, and analysis name. See createParams.m for full
%   list.

%% Pre-Analysis Setup

% Add CoSMoMVPA to the MATLAB search path
addpath(genpath('/path/to/CoSMoToolbox'));

% Required for boostrap code to function
%addpath('/path/to/CANLab-Multivariate-Scripts/functions/cosmo_crossvalidate_bootstrap.m');
addpath(fileparts(which('cosmo_crossvalidate_bootstrap.m')));

% turn cosmo warnings off
%cosmo_warning('off');

%% Set Analysis Parameters & Paths
% Load subject IDs, ROIs, and Condition flags
if exist('flag','var')==0
    
    %Select parameter file is flag does not exist
    uiopen('*.mat')
    
end

% Filepath for results folder
parentDir = fileparts(study_path);

% Base output directory name
switch analysisType
    case 'Searchlight'
        analysis=strcat(parentDir,'/',analysisName,'_',classType,'_',...
            analysisType,'_',trialAnalysis,'_',metric,'_',...
            num2str(searchlightSize),'_',conds{1,1},'_',conds{1,2});
    otherwise
        if exist('subConds','var')
            analysis=strcat(parentDir,'/',analysisName,'_',classType,'_',...
                analysisType,'_',trialAnalysis,'_',conds{1,1},'_',...
                conds{1,2},'_',subConds{1,1},'_',subConds{1,2});
        else
            analysis=strcat(parentDir,'/',analysisName,'_',classType,'_',...
                analysisType,'_',trialAnalysis,'_',conds{1,1},'_',conds{1,2});
        end
end

% Bootstrap flag
switch bootstrap.flag
    case 'Yes'
        analysis = [analysis '_Bootstrap'];
end

%% Main Body
for iteration=1:length(subjects)*length(rois)
    
    % Loop count for all subject/region combinations
    if iteration==-1
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
    output_path = fullfile(out_path, subject);
    spm_fn = [data_path '/SPM.mat'];
    
    % create the output path if it doesn't already exist
    if ~exist(output_path, 'dir')
        mkdir(output_path)
    end
    
    % Path to current region mask
    curROI = fullfile(roi_path, ROI);
    
    % Current region name
    regionName=erase(ROI,{'reslice_','_bilat.nii','-'});
    
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
    
    %%% Test/Train Flag
    % If set to 'Trial', chunks will be reset by trial instead of by run.
    % Each trials is now its own separate chunk
    switch trialAnalysis
        case 'Trial'
            currDataset.sa.chunks=(1:1:length(currDataset.sa.chunks))';
    end
    
    %%% Indentify trials of interest for each condition
    if exist('subConds','var')
        
        for ii=1:length(conds)
            
            subCond.(conds{1,ii})=contains(currDataset.sa.labels, conds{1,ii});
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
    fprintf('Number of samples: %i\n',size(currDataset.samples,1));
    fprintf('Number of features (voxels): %i\n',size(currDataset.samples,2));
    fprintf('Number of chunks (runs): %i\n',length(unique(currDataset.sa.chunks)));
    
    %% classifier ROI analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define which classifier to use, using a function handle.
    % Alternatives are @cosmo_classify_{svm,matlabsvm,libsvm,nn,naive_bayes}
    
    % Set partition scheme. odd_even is fast; for publication-quality analysis
    % nfold_partitioner (cross-validation) is recommended.
    % Alternatives are:
    % - cosmo_nfold_partitioner    (take-one-chunk-out crossvalidation)
    % - cosmo_nchoosek_partitioner (take-K-chunks-out  "             ").
    % We will also make sure the partitions are *balanced* (targets are evenly
    % distributed).
    partitions = cosmo_nfold_partitioner(currDataset);
    
    %%% Test/Train Flag
    % If set to 'Trial', partitions will not be balanced. Balancing is not
    % possible for 'Trial' as n-1 trials are used to train and tested on
    % the remaining trial (e.g. 100 Trials = train on 99, test on 1)
    switch trialAnalysis
        case 'Run'
            partitions = cosmo_balance_partitions(partitions, currDataset);
    end
    
    fprintf('There are %d partitions\n', numel(partitions.train_indices));
    fprintf('# train samples:%s\n', sprintf(' %d', cellfun(@numel, ...
        partitions.train_indices)));
    fprintf('# test samples:%s\n', sprintf(' %d', cellfun(@numel, ...
        partitions.test_indices)));
    
    % Set any other options - see help cosmo_crossvalidate
    opt = struct();
    
    switch analysisType
        case 'Searchlight'
            % This analysis identifies brain regions that discriminate the categories
            % using a whole-brain searchlight procedure. The classifier uses an n-fold
            % partitioning scheme (cross-validation) and a Linear Discriminant Analysis
            % (LDA) classifier.
            
            % Use the cosmo_cross_validation_measure and set its parameters
            % (classifier and partitions) in a measure_args struct.
            measure = @cosmo_crossvalidation_measure;
            
            % Defines which classifier to use, using a function handle.
            classifier.function = @cosmo_classify_lda;
            classifier.name = 'LDA';
            opt.classifier = classifier.function;
            
            % Set partition scheme.
            opt.partitions = partitions;
            
            try
                % Define a neighborhood with approximately 100 voxels in 
                % each searchlight.
                nbrhood=cosmo_spherical_neighborhood(currDataset,metric,...
                    searchlightSize);
                
                %%% Test/Train Flag
                % If set to 'Trial', partitions are not checked for balance
                % because they cannot be. Code will crash/not classify if
                % balancing is not set to be ignored.
                switch trialAnalysis
                    case 'Trial'
                        opt.check_partitions=false;
                end
                
                % Run the searchlight
                searchResults = cosmo_searchlight...
                    (currDataset,nbrhood,measure,opt);
                save([output_path '/searchlightResults_' regionName '_' metric '_' ...
                    num2str(searchlightSize) '.mat'],'searchResults');
                
                % print output dataset
                fprintf('Dataset output:\n');
                
                % Define output location
                if ~exist('subConds','var')
                    output=strcat(output_path,'/',subject,'_',...
                        classifier.name,'_Searchlight_',regionName,'_',...
                        metric,'_',num2str(searchlightSize),'_',...
                        conds{1,1},'_vs_',conds{1,2},'.nii');
                else
                    output=strcat(output_path,'/',subject,'_',...
                        classifier.name,'_Searchlight_',regionName,'_',...
                        metric,'_',num2str(searchlightSize),'_',...
                        conds{1,1},'_vs_',conds{1,2},'_',subConds{1,1},'.nii');
                end
                
                % Store results to disc
                cosmo_map2fmri(searchResults, output);
                
                % Gunzip output to save space
                setenv('output',output);
                !gzip $output
            catch ME
                disp([ME.message ' Unable to create searchlight for ' ...
                    regionName '!']);
            end
            
        otherwise
            % This analysis calculates classification accuracy using all
            % voxels within the ROI. The classifier uses an n-fold
            % partitioning scheme (cross-validation) and a
            % Support Vector Machine (SVM) classifier.
            
            % Defines which classifier to use, using a function handle.
            classifier.function = @cosmo_classify_svm;
            classifier.name = 'SVM';
            
            % Set any other options - see help cosmo_crossvalidate
            opt.normalization = 'zscore';
            opt.max_feature_count = size(currDataset.samples,2);
            
            %%% Test/Train Flag
            % If set to 'Trial', partitions are not checked for balance
            % because they cannot be. Code will crash/not classify if
            % balancing is not set to be ignored.
            switch trialAnalysis
                case 'Trial'
                    opt.check_partitions=false;
            end
            
            %Bootstrapping
            switch bootstrap.flag
                case 'Yes'
                    opt.check_partitions=false;
                    
                    for i=1:bootstrap.numRuns
                        
                        runTrials=find(currDataset.sa.chunks==i);
                        
                        if i==1
                            diff = length(runTrials)-...
                                length(partitions.test_indices{i});
                            counter=0;
                        end
                        
                        if ~isempty(runTrials)
                            block(i).begin = min(runTrials)-counter;
                            counter=counter+diff;
                            block(i).end = max(runTrials)-counter;
                        end
                        
                    end
                    
                    if iteration==1
                        try
                            permFiles=dir([parentDir '/*permutation*.mat']);
                            for i = 1:length(permFiles)
                                dates(i) = datenum(permFiles(i).date);
                            end
                            
                            [tmp,i] = max(dates);
                            load([parentDir filesep permFiles(i).name]);
                            
                        catch
                            
                            for i=1:bootstrap.perm
                                for j=1:length(block)
                                    permutation(i).(...
                                        bootstrap.structNames{j}) = ...
                                        randperm(length(bootstrap.trialsPerRun));
                                end
                            end
                            
                            %Save Date/Time in filename
                            timePermute = datestr(now,'yyyymmddHHMMSS');
                            save([out_path '_' classType '_permutation_'...
                                timePermute '.mat'],'permutation');
                        end
                    end
                    
                    for i=1:bootstrap.perm
                        
                        if i==1
                            regionName=erase(ROI,{'reslice_','_bilat.nii'});
                        end
                        currPermute = permutation(i);
                        
                        try
                            [predictions, accuracy] = ...
                                cosmo_crossvalidate_bootstrap...
                                (currDataset,classifier.function,...
                                partitions,opt,bootstrap.structNames,...
                                currPermute,block);
                        catch
                            accuracy = NaN;
                        end
                        finalPred.(regionName)(:,i)=predictions;
                        finalAcc(i,1)=accuracy;
                        
                    end
                    concatAccuracy(:,regionCount) = finalAcc;
                    
                    % Run Standard SVM classification with no bootstrap
            end
            
            %regionName=erase(ROI,{'reslice_','_bilat.nii','-'});
            try
                [predictions,accuracy] = cosmo_crossvalidate...
                    (currDataset,classifier.function,partitions,opt);
                % Report results in command window
                fprintf('Accuracy: %0.2f\n',accuracy);
            catch
                fprintf('Classification Failed!\n');
                predictions(1:length(currDataset.sa.chunks),1) = NaN;
                accuracy = NaN;
            end
            
            finalPredictions(:,regionCount) = predictions;
            clear predictions;
    end
    
    %% Calculating Classification Accuracy
    if strcmpi(analysisType,'Searchlight')==0
        
        % Caluclate accuracy by individual trial
        switch bootstrap.flag
            case 'Yes'
                for j=1:size(finalPred.(regionName),2)
                    for k=1:size(finalPred.(regionName),1)
                        if isnan(finalPred.(regionName)(k,j))==1
                            correct.(regionName)(k,j)=NaN;
                        elseif finalPred.(regionName)(k,j)==currDataset.sa.targets(k)
                            correct.(regionName)(k,j)=1;
                        else
                            correct.(regionName)(k,j)=0;
                        end
                    end
                end
            case 'No'
                for j=1:length(finalPredictions(:,regionCount))
                    if isnan(finalPredictions(j,regionCount))==1
                        correct.(regionName)(j,1)=NaN;
                    elseif finalPredictions(j,regionCount)==currDataset.sa.targets(j)
                        correct.(regionName)(j,1)=1;
                    else
                        correct.(regionName)(j,1)=0;
                    end
                end
        end
        
        % Calculate accuracy by behavioral trial type (e.g., Hit acc, Miss acc)
        
        % Pull trial information from ROI by target condition
        Cond1 = ~cellfun(@isempty, strfind(currDataset.sa.labels, conds{1,1}));
        Cond1idx = find(Cond1 == 1);
        Cond2 = ~cellfun(@isempty, strfind(currDataset.sa.labels, conds{1,2}));
        Cond2idx = find(Cond2 == 1);
        trialInfo.Cond1 = currDataset.sa.labels(Cond1idx);
        trialInfo.Cond2 = currDataset.sa.labels(Cond2idx);
        
        switch bootstrap.flag
            case 'No'
                
                % Pull accuracy information by target condition
                acc.Cond1 = correct.(regionName)(Cond1idx');
                acc.Cond2 = correct.(regionName)(Cond2idx');
                acc.Cond1 = acc.Cond1';
                acc.Cond2 = acc.Cond2';
                
                % Find behavioral trials for each condition
                acc.rec1 = ~cellfun(@isempty, strfind(trialInfo.Cond1,'Hit'));
                acc.fam1 = ~cellfun(@isempty, strfind(trialInfo.Cond1,'Familiar'));
                acc.miss1 = ~cellfun(@isempty, strfind(trialInfo.Cond1,'Miss'));
                
                acc.rec2 = ~cellfun(@isempty, strfind(trialInfo.Cond2,'Hit'));
                acc.fam2 = ~cellfun(@isempty, strfind(trialInfo.Cond2,'Familiar'));
                acc.miss2 = ~cellfun(@isempty, strfind(trialInfo.Cond2,'Miss'));
                
                acc.recFinal1 = acc.Cond1(acc.rec1==1);
                acc.famFinal1 = acc.Cond1(acc.fam1==1);
                acc.missFinal1 = acc.Cond1(acc.miss1==1);
                
                acc.recFinal2 = acc.Cond2(acc.rec2==1);
                acc.famFinal2 = acc.Cond2(acc.fam2==1);
                acc.missFinal2 = acc.Cond2(acc.miss2==1);
                
                if iteration==1
                    finalTableTrialTypeAcc=cell(length(subjects)+1,length(rois)+2);
                    finalTableTrialTypeAcc{1,1}='subjectid';
                    finalTableTrialTypeAcc{1,2}='Trial Type';
                    tempcount=3;
                    
                    for headerTrialType=1:length(rois)
                        finalTableTrialTypeAcc{1,tempcount}=strcat...
                            (rois{1,headerTrialType}(1:end-4),'_trainedRecAccuracy');
                        finalTableTrialTypeAcc{1,tempcount+1}=strcat...
                            (rois{1,headerTrialType}(1:end-4),'_novelRecAccuracy');
                        
                        finalTableTrialTypeAcc{1,tempcount+2}=strcat...
                            (rois{1,headerTrialType}(1:end-4),'_trainedRecNumTrials');
                        finalTableTrialTypeAcc{1,tempcount+3}=strcat...
                            (rois{1,headerTrialType}(1:end-4),'_novelRecNumTrials');
                        
                        finalTableTrialTypeAcc{1,tempcount+4}=strcat...
                            (rois{1,headerTrialType}(1:end-4),'_trainedFamAccuracy');
                        finalTableTrialTypeAcc{1,tempcount+5}=strcat...
                            (rois{1,headerTrialType}(1:end-4),'_novelFamAccuracy');
                        
                        finalTableTrialTypeAcc{1,tempcount+6}=strcat...
                            (rois{1,headerTrialType}(1:end-4),'_trainedFamNumTrials');
                        finalTableTrialTypeAcc{1,tempcount+7}=strcat...
                            (rois{1,headerTrialType}(1:end-4),'_novelFamNumTrials');
                        
                        finalTableTrialTypeAcc{1,tempcount+8}=strcat...
                            (rois{1,headerTrialType}(1:end-4),'_trainedMissAccuracy');
                        finalTableTrialTypeAcc{1,tempcount+9}=strcat...
                            (rois{1,headerTrialType}(1:end-4),'_novelMissAccuracy');
                        
                        finalTableTrialTypeAcc{1,tempcount+10}=strcat...
                            (rois{1,headerTrialType}(1:end-4),'_trainedMissNumTrials');
                        finalTableTrialTypeAcc{1,tempcount+11}=strcat...
                            (rois{1,headerTrialType}(1:end-4),'_novelMissNumTrials');
                        
                        tempcount=tempcount+12;
                    end
                    
                    row=2;
                    headerTrialType=3;
                    clear tempcount;
                    TrialTypeCombo = {strcat(conds{1,1},'_v_',conds{1,2})};
                end
                
                % Add subject, trial type, and accuracy to table
                finalTableTrialTypeAcc{row,1}=subject;
                finalTableTrialTypeAcc{row,2}=TrialTypeCombo{1,1};
                finalTableTrialTypeAcc{row,headerTrialType}=num2str...
                    (nanmean(acc.recFinal1));
                headerTrialType=headerTrialType+1;
                finalTableTrialTypeAcc{row,headerTrialType}=num2str...
                    (nanmean(acc.recFinal2));
                headerTrialType=headerTrialType+1;
                
                finalTableTrialTypeAcc{row,headerTrialType}=num2str...
                    (length(acc.recFinal1) - sum(isnan(acc.recFinal1)));
                headerTrialType=headerTrialType+1;
                finalTableTrialTypeAcc{row,headerTrialType}=num2str...
                    (length(acc.recFinal2) - sum(isnan(acc.recFinal2)));
                headerTrialType=headerTrialType+1;
                
                finalTableTrialTypeAcc{row,headerTrialType}=num2str...
                    (nanmean(acc.famFinal1));
                headerTrialType=headerTrialType+1;
                finalTableTrialTypeAcc{row,headerTrialType}=num2str...
                    (nanmean(acc.famFinal2));
                headerTrialType=headerTrialType+1;
                
                finalTableTrialTypeAcc{row,headerTrialType}=num2str...
                    (length(acc.famFinal1) - sum(isnan(acc.famFinal1)));
                headerTrialType=headerTrialType+1;
                finalTableTrialTypeAcc{row,headerTrialType}=num2str...
                    (length(acc.famFinal2) - sum(isnan(acc.famFinal2)));
                headerTrialType=headerTrialType+1;
                
                finalTableTrialTypeAcc{row,headerTrialType}=num2str...
                    (nanmean(acc.missFinal1));
                headerTrialType=headerTrialType+1;
                finalTableTrialTypeAcc{row,headerTrialType}=num2str...
                    (nanmean(acc.missFinal2));
                headerTrialType=headerTrialType+1;
                
                finalTableTrialTypeAcc{row,headerTrialType}=num2str...
                    (length(acc.missFinal1) - sum(isnan(acc.missFinal1)));
                headerTrialType=headerTrialType+1;
                finalTableTrialTypeAcc{row,headerTrialType}=num2str...
                    (length(acc.missFinal2) - sum(isnan(acc.missFinal2)));
                headerTrialType=headerTrialType+1;
                
        end
        
    end
    
    %% Save text output of SVM Classification
    try
        if strcmpi(analysisType,'Searchlight')==0
            
            % Create a tidyverse formatted table for final statistical analysis
            TrialTypeCombo = {strcat(conds{1,1},'_v_',conds{1,2})};
            
            % create subjectid and roiid columns
            subjectid   = repmat({subject}, length(TrialTypeCombo), 1);
            roiid       = repmat({regionName}, length(TrialTypeCombo), 1);
            
            % create the stats table
            stats_table = table(subjectid, roiid, TrialTypeCombo, accuracy);
            
            % write the stats table
            filename = sprintf('sub-%s_roi-%s_%s_statistics-table.csv',...
                subject, regionName, classifier.name);
            writetable(stats_table, fullfile(output_path, filename));
            
            %% Create aggregate table of classification
            % Create headers
            if iteration==1
                finalTable=cell(length(subjects)+1,length(rois)+2);
                finalTable{1,1}='subjectid';
                finalTable{1,2}='Trial Type';
                tempcount=3;
                switch bootstrap.flag
                    case 'No'
                        for header=1:length(rois)
                            finalTable{1,tempcount}=strcat...
                                (rois{1,header}(1:end-4),'_Accuracy');
                            tempcount=tempcount+1;
                        end
                    case 'Yes'
                        for header=1:length(rois)
                            finalTable{1,tempcount}=strcat...
                                (rois{1,header}(1:end-4),'_true_accuracy');
                            finalTable{1,tempcount+1}=strcat...
                                (rois{1,header}(1:end-4),'_perm_skew');
                            finalTable{1,tempcount+2}=strcat...
                                (rois{1,header}(1:end-4),'_perm_kurtosis');
                            tempcount=tempcount+3;
                        end
                end
                
                finalTable{1,tempcount}='numTargets';
                finalTable{1,tempcount+1}='numTrain';
                
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
            switch bootstrap.flag
                case 'No'
                    finalTable{row,header}=accuracy;
                    header=header+1;
                case 'Yes'
                    save([output_path filesep ROI '_permutedAcc.mat'],'finalAcc');
                    finalTable{row,header}=accuracy;
                    finalTable{row,header+1}=num2str(skewness(finalAcc));
                    finalTable{row,header+2}=num2str(kurtosis(finalAcc));
                    header=header+3;
                    
                    clear finalAcc;
            end
            
            % Drops to next row if remainder is 0 (e.g. all regions have been
            % entered for a given subject)
            if iterCheck == 0
                finalTable{row,header}=num2str(...
                    length(partitions.test_indices{1,1}));
                finalTable{row,header+1}=num2str(...
                    length(partitions.train_indices{1,1}));
                row=row+1;
                header=3;
                headerTrialType=3;
                
                switch bootstrap.flag
                    case 'No'
                        % Save mat file with predictions separately
                        save([output_path filesep 'finalPredictions.mat'],...
                            'finalPredictions');
                        save([output_path filesep 'finalAccuracy.mat'],'correct');
                        clear finalPredictions correct;
                    case 'Yes'
                        save([output_path filesep 'allROIPermuteACC.mat'],...
                            'concatAccuracy');
                        save([output_path filesep 'finalPredAllROIs.mat'],...
                            'finalPred');
                        clear concatAccuracy finalPred;
                end
            end
        end
    catch
        error(['Unable to record classifier accuracy for ' subject '!']);
    end
    
end

%% Save summary files of SVM classifation.
try
    if strcmpi(analysisType,'Searchlight')==0
        
        % Save classifier outputs
        save([fileparts(output_path) filesep 'finalTable.mat'],'finalTable');
        switch bootstrap.flag
            case 'No'
                save([fileparts(output_path) filesep 'finalTableTrialType.mat'],...
                    'finalTableTrialTypeAcc')
        end
        
        % Copy summary file to temp variable and remove headers
        temp=finalTable;
        temp(:,1:2)=[];
        temp(1,:)=[];
        
        % Remove target/trial number columns off end
        temp(:,end-1:end)=[];
        
        % Create table of double type values for statistical tests
        tempDouble=zeros(size(temp));
        
        for i=1:size(temp,1)
            for ii=1:size(temp,2)
                tempDouble(i,ii)=double(temp{i,ii});
            end
        end
        
        clear temp;
        
        switch bootstrap.flag
            case 'No'
                %% One-Sample T-test
                % Subject Counter
                index=size(finalTable,1);
                
                for i=2:(size(finalTable,2)-2)
                    % Row header information
                    if i==2
                        
                        labelVec={strcat('Statistics (Tested against chance: ',num2str...
                            (1/length(conds)),')'),'Significant at p<0.05?','mean',...
                            'min','max','p value','Confidence Interval','t value',...
                            'df','sd','se'};
                        
                        finalTable(index+1:index+length(labelVec),1)=labelVec;
                        
                    else
                        % ROI specific statistics
                        [h,p,ci,stats]=ttest(tempDouble(:,i-2),(1/length(conds)));
                        if h==1
                            sig='Yes';
                        else
                            sig='No';
                        end
                        
                        subjVec={sig,mean(tempDouble(:,i-2)),min(tempDouble(:,i-2)),...
                            max(tempDouble(:,i-2)),p,ci,stats.tstat,num2str(stats.df),...
                            stats.sd,stats.sd/sqrt(size(tempDouble,1))};
                        
                        finalTable(index+2:index+length(subjVec)+1,i)=subjVec;
                        
                        clear h p ci stats sig;
                    end
                    
                    % Save mat file with statistics separately
                    save([out_path '_' classType filesep 'finalTable_stats.mat'],'finalTable');
                    
                end
            case 'Yes'
                %% Compare True Accuracy to Permutation Classification
                criticalCutoff = mean(tempDouble);
                
                for i=1:length(rois)
                    
                    % Create aggregate table of classification
                    if i==1
                        finalTableBootstrap=cell(4,length(rois)+1);
                        finalTableBootstrap{2,1}='True Accuracy';
                        finalTableBootstrap{3,1}='Permutations greater than True Accuracy';
                        finalTableBootstrap{4,1}='Bootstrap p value';
                        tempcount=2;
                        
                        for header=1:length(rois)
                            finalTableBootstrap{1,tempcount}=strcat...
                                (rois{1,header}(1:end-4),'_Accuracy');
                            tempcount=tempcount+1;
                        end
                        
                        row=2;
                        header=3;
                        clear tempcount;
                    end
                    
                    % Aggregate final accuracies for all permuatations for all subjects
                    for j=1:length(subjects)
                        
                        file = dir([fileparts(output_path) '/' subjects{j} '/' rois{i} '*']);
                        load([file.folder '/' file.name]);
                        
                        tempROI(:,j) = finalAcc;
                        
                        clear finalAcc;
                        
                    end
                    
                    % Average all accuracies per permutation and sort descending
                    avgAcc = mean(tempROI,2);
                    avgAcc = sort(avgAcc,'descend');
                    
                    % Flag any accuracy great than or equal to the true accuracy
                    counter = 0;
                    for j=1:length(avgAcc)
                        if avgAcc(j) >= criticalCutoff(i)
                            counter=counter+1;
                        end
                    end
                    
                    % Add values to final table
                    finalTableBootstrap{2,i+1} = num2str(criticalCutoff(i));
                    finalTableBootstrap{3,i+1} = num2str(counter);
                    finalTableBootstrap{4,i+1} = num2str(counter/length(avgAcc));
                    
                    clear avgAcc counter;
                    
                end
        end
        
        
        %% Write output summary file
        
        switch bootstrap.flag
            case 'No'
                % Save classifier accuracy and one-sample t-tets results to
                % CSV file
                file = fopen([out_path '_MVPA' filesep  'allAccuraciesSummary.csv'], 'w');
                
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
                
                % Save classifier accuracy by trial type/behavior to CSV
                file = fopen([out_path '_MVPA' filesep 'allAccuraciesByTrialType.csv'], 'w');
                
                for a=1:size(finalTableTrialTypeAcc,1)
                    for b=1:size(finalTableTrialTypeAcc,2)
                        var = eval('finalTableTrialTypeAcc{a,b}');
                        try
                            fprintf(file, '%s', var);
                        end
                        fprintf(file, ',');
                    end
                    fprintf(file, '\n');
                end
                fclose(file);
                
                %Run R Markdown graph
                try
                    % Directory variables
                    curDir = pwd;
                    markdown = which('anovaGraphOnlyMarkdown.r');
                    
                    % Call R script from command line
                    cmd = ['Rscript -e "library(knitr); dataPath <- ''' analysis...
                        '''; knit(''' markdown ''', ''anovaGraphOnlyMarkdown.html'')"'];
                    [status,cmdout] = system(cmd);
                    
                    % Move HTML to output folder
                    setenv('dataPath',analysis);
                    setenv('curDir',curDir);
                    !mv -t $dataPath $curDir/analyses/anovaGraphOnlyMarkdown.html figure/
                    
                catch
                    message('Unable to save table and figure!');
                end
                
            case 'Yes'
                % Save comparison of permuted and true accuracy to CSV file
                file = fopen([fileparts(output_path) filesep  'permAccuraciesSummary.csv'], 'w');
                
                for a=1:size(finalTableBootstrap,1)
                    for b=1:size(finalTableBootstrap,2)
                        var = eval('finalTableBootstrap{a,b}');
                        try
                            fprintf(file, '%s', var);
                        end
                        fprintf(file, ',');
                    end
                    fprintf(file, '\n');
                end
                fclose(file);
        end
        clear;
    end
end
