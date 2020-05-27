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

%% Set Analysis Parameters & Paths
% Load all relevent project information
if exist('flag','var') == 0
    
    %Select parameter file is flag does not exist
    [file,path]=uigetfile('*.mat','Select params file');
    filename=fullfile(path,file);
    load(filename);
    
end

% Required for boostrap code to function
addpath(fileparts(which('cosmo_crossvalidate_bootstrap.m')));

% turn cosmo warnings off
cosmo_warning('off');

% Filepath for results folder
parentDir = directory.Model;

% Base output directory name
analysis = [directory.Analysis filesep 'models' filesep file(1:end-4)];

%Debug
subjects(2)=[];

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
    masks=masks(9:13);
    bootstrap.perm = 500;
    
    for curMask = 1:length(masks)
        
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
                files = dir([studyPath filesep subject filesep 'Run*']);
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
        
        %%% Test/Train Flag %%%
        % If set to 'Trial', chunks will be reset by trial instead of by run.
        % Each trials is now its own separate chunk
        %
        % This has been removed from the parameters but can be included by
        % uncommenting this code and resetting the 'trialAnalysis'
        % variable to 'Trial'.
        %         switch trialAnalysis
        %             case 'Run'
%         switch trialAnalysis
%             case 'Trial'
%                 currDataset.sa.chunks=(1:1:length(currDataset.sa.chunks))';
%         end
        
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
        %
        % This has been removed from the parameters but can be included by
        % uncommenting this code and resetting the 'trialAnalysis'
        % variable to anything.
        %         switch trialAnalysis
        %             case 'Run'
        partitions = cosmo_balance_partitions(partitions, currDataset);
        %         end
        
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
                    save([outputPath '/searchlightResults_' regionName '_' metric '_' ...
                        num2str(searchlightSize) '.mat'],'searchResults');
                    
                    % print output dataset
                    fprintf('Dataset output:\n');
                    
                    % Define output location
                    if ~exist('subConds','var')
                        output=strcat(outputPath,'/',subject,'_',...
                            classifier.name,'_Searchlight_',regionName,'_',...
                            metric,'_',num2str(searchlightSize),'_',...
                            taskInfo.Conditions{1},'_vs_',taskInfo.Conditions{2},'.nii');
                    else
                        output=strcat(outputPath,'/',subject,'_',...
                            classifier.name,'_Searchlight_',regionName,'_',...
                            metric,'_',num2str(searchlightSize),'_',...
                            taskInfo.Conditions{1},'_vs_',taskInfo.Conditions{2},'_',subConds{1,1},'.nii');
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
                %
                % This has been removed from the parameters but can be 
                % included by uncommenting this code and resetting the 
                % 'trialAnalysis' variable to 'Trial.
%                 switch trialAnalysis
%                     case 'Trial'
%                         opt.check_partitions=false;
%                 end
                
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
                                permFiles=dir([fileparts(parentDir) '/*permutation*.mat']);
                                for i = 1:length(permFiles)
                                    dates(i) = datenum(permFiles(i).date);
                                end
                                
                                [tmp,i] = max(dates);
                                load([fileparts(parentDir) filesep permFiles(i).name]);
                                
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
                                save([analysis '_' classType '_permutation_'...
                                    timePermute '.mat'],'permutation');
                            end
                        end
                        
                        for i=1:bootstrap.perm
                            
                            currPermute = permutation(i);
                            
                            % Unequal Trial Check  %%% IN PROGRESS %%%
%                             for j=1:taskInfo.Runs
%                                 if length(partitions.test_indices{j}) ~= length(permutation(i).(['perm' num2str(j)]))
%                                     [~,I] = max(permutation(i).(['perm' num2str(j)]));
%                                     permutation(i).(['perm' num2str(j)])(I)=[];
%                                     [~,I] = max(permutation(i).(['perm' num2str(j)]));
%                                     permutation(i).(['perm' num2str(j)])(I)=[];
%                                 end
%                             end
                            
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
                        concatAccuracy(:,curMask) = finalAcc;
                        
                end
                
                % Run Standard SVM classification with no bootstrap
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
                
                finalPredictions(:,curMask) = predictions;
                clear predictions;
        end
        
        %% Calculating Classification Accuracy
        if strcmpi(analysisType,'Searchlight')==0
            
            % Caluclate classifier accuracy by individual trial
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
                    for j=1:length(finalPredictions(:,curMask))
                        if isnan(finalPredictions(j,curMask))==1
                            correct.(regionName)(j,1)=NaN;
                        elseif finalPredictions(j,curMask)==currDataset.sa.targets(j)
                            correct.(regionName)(j,1)=1;
                        else
                            correct.(regionName)(j,1)=0;
                        end
                    end
            end
            
        end
        
        %% Save text output of SVM Classification
        try
            if strcmpi(analysisType,'Searchlight')==0
                
                % Create a tidyverse formatted table for final statistical analysis
                for i=1:length(taskInfo.Conditions)
                    if i==1
                        TrialTypeCombo = taskInfo.Conditions{i};
                    else
                        TrialTypeCombo = strcat(TrialTypeCombo,'_v_',taskInfo.Conditions{i});
                    end
                end
                
                TrialTypeCombo = {TrialTypeCombo};
                
                % create subjectid and roiid columns
                subjectid   = repmat({subjects{iteration}}, ...
                    length(TrialTypeCombo), 1);
                roiid       = repmat({regionName}, ...
                    length(TrialTypeCombo), 1);
                
                % create the stats table
                stats_table = table(subjectid, roiid, TrialTypeCombo, accuracy);
                
                % write the stats table
                filename = sprintf('sub-%s_roi-%s_%s_statistics-table.csv',...
                    subjects{iteration}, regionName, classifier.name);
                writetable(stats_table, fullfile(outputPath, filename));
                
                %% Create aggregate table of classification
                % Create headers
                if iteration==1 && curMask==1
                    summary=cell(length(subjects)+1,length(masks)+2);
                    summary{1,1}='subjectid';
                    summary{1,2}='Trial Type';
                    tmpCnt=3;
                    switch bootstrap.flag
                        case 'No'
                            for header=1:length(masks)
                                summary{1,tmpCnt}=[masks(header).name(1:end-7)...
                                    '_classAcc'];
                                tmpCnt=tmpCnt+1;
                            end
                        case 'Yes'
                            for header=1:length(masks)
                                summary{1,tmpCnt}=[masks(header).name(1:end-7)...
                                    '_true_accuracy'];
                                summary{1,tmpCnt+1}=[masks(header).name(1:end-7)...
                                    '_perm_skew'];
                                summary{1,tmpCnt+2}=[masks(header).name(1:end-7)...
                                    '_perm_kurtosis'];
                                tmpCnt=tmpCnt+3;
                            end
                    end
                    
                    summary{1,tmpCnt}='numTargets';
                    summary{1,tmpCnt+1}='numTrain';
                    
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
                switch bootstrap.flag
                    case 'No'
                        summary{row,header}=accuracy;
                        header=header+1;
                    case 'Yes'
                        save([outputPath filesep regionName...
                            '_permutedAcc.mat'],'finalAcc');
                        summary{row,header}=accuracy;
                        summary{row,header+1}=num2str(skewness(finalAcc));
                        summary{row,header+2}=num2str(kurtosis(finalAcc));
                        header=header+3;
                        
                        clear finalAcc;
                end
                
                % Drops to next row if remainder is 0 (e.g. all regions have been
                % entered for a given subject)
                if iterCheck == 0
                    summary{row,header}=num2str(...
                        length(partitions.test_indices{1,1}));
                    summary{row,header+1}=num2str(...
                        length(partitions.train_indices{1,1}));
                    row=row+1;
                    header=3;
                    headerTrialType=3;
                    
                    switch bootstrap.flag
                        case 'No'
                            % Save mat file with predictions separately
                            save([outputPath filesep 'finalPredictions.mat'],...
                                'finalPredictions');
                            save([outputPath filesep 'finalAccuracy.mat'],'correct');
                            clear finalPredictions correct;
                        case 'Yes'
                            save([outputPath filesep 'allROIPermuteACC.mat'],...
                                'concatAccuracy');
                            save([outputPath filesep 'finalPredAllROIs.mat'],...
                                'finalPred');
                            clear concatAccuracy finalPred;
                    end
                end
            end
        catch
            error(['Unable to record classifier accuracy for ' subjects{iteration} '!']);
        end
        
    end
end

%% Save summary files of SVM classifation.
try
    if strcmpi(analysisType,'Searchlight')==0
        
        % Save summary classifier accuracies as mat file
        save([fileparts(outputPath) filesep 'summary.mat'],'summary');
        
        % Copy summary file to temp variable and remove headers
        temp=summary;
        temp(:,1:2)=[];
        temp(1,:)=[];
        
        % Remove target/trial number columns off end
        temp(:,end-1:end)=[];
        
        % Create table of double type values for statistical tests
        tempDouble=zeros(size(temp));
        
        for i=1:size(temp,1)
            for ii=1:size(temp,2)
                if isnumeric(temp{i,ii})==0
                    %tempDouble(i,ii)=double(temp{i,ii});
                    tempDouble(i,ii)=str2double(temp{i,ii});
                else
                    tempDouble(i,ii)=temp{i,ii};
                end
            end
        end
        
        clear temp;
        
        switch bootstrap.flag
            case 'No'
                %% One-Sample T-test
                % Subject Counter
                index=size(summary,1);
                
                for i=2:(size(summary,2)-2)
                    % Row header information
                    if i==2
                        
                        labelVec={strcat('Statistics (Tested against chance: ',num2str...
                            (1/length(taskInfo.Conditions)),')'),...
                            'Significant at p<0.05?','mean','min','max',...
                            'p value','Confidence Interval','t value',...
                            'df','sd','se'};
                        
                        summary(index+1:index+length(labelVec),1)=labelVec;
                        
                    else
                        % ROI specific statistics
                        [h,p,ci,stats]=ttest(tempDouble(:,i-2),...
                            (1/length(taskInfo.Conditions)));
                        if h==1
                            sig='Yes';
                        else
                            sig='No';
                        end
                        
                        subjVec={sig,mean(tempDouble(:,i-2)),min(tempDouble(:,i-2)),...
                            max(tempDouble(:,i-2)),p,ci,stats.tstat,num2str(stats.df),...
                            stats.sd,stats.sd/sqrt(size(tempDouble,1))};
                        
                        summary(index+2:index+length(subjVec)+1,i)=subjVec;
                        
                        clear h p ci stats sig;
                    end
                    
                    % Save mat file with statistics separately
                    save([fileparts(outputPath) filesep 'summary_stats.mat'],'summary');
                    
                end
            case 'Yes'
                %% Compare True Accuracy to Permutation Classification
                criticalCutoff = mean(tempDouble);
                
                for i=1:length(subjects)
                    files = dir([fileparts(outputPath) '/' subjects{i} '/*permutedAcc.mat']);
                    
                    for j=1:length(files)
                        
                        % Create aggregate table of classification
                        if i==1 && j==1
                            finalTableBootstrap=cell(4,length(files)+1);
                            finalTableBootstrap{2,1}='True Accuracy';
                            finalTableBootstrap{3,1}='Permutations greater than True Accuracy';
                            finalTableBootstrap{4,1}='Bootstrap p value';
                            tmpCnt=2;
                        
                            for header=1:length(masks)
                                finalTableBootstrap{1,tmpCnt}=strcat...
                                    (files(header).name(1:end-7),'_classAcc');
                                tmpCnt=tmpCnt+1;
                            end
                        
                            row=2;
                            header=3;
                            clear tempcount;
                            
                        end
                        
                        % Aggregate final accuracies for all permuatations for all subjects
                        load([files(j).folder '/' files(j).name]);
                        tempROI(:,j) = finalAcc;
                        clear finalAcc;
                        
                    end
                    
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
                % Save classifier accuracy and one-sample t-test results to
                % CSV file
                file = fopen([fileparts(outputPath) filesep  'allClassAccSummary.csv'], 'w');
                
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
                
                %Run R Markdown graph %%% IN DEVELOPMENT %%%
%                 try
%                     % Directory variables
%                     curDir = pwd;
%                     markdown = which('anovaGraphOnlyMarkdown.r');
%                     
%                     % Call R script from command line
%                     cmd = ['Rscript -e "library(knitr); dataPath <- ''' analysis...
%                         '''; knit(''' markdown ''', ''anovaGraphOnlyMarkdown.html'')"'];
%                     [status,cmdout] = system(cmd);
%                     
%                     % Move HTML to output folder
%                     setenv('dataPath',analysis);
%                     setenv('curDir',curDir);
%                     !mv -t $dataPath $curDir/analyses/anovaGraphOnlyMarkdown.html figure/
%                     
%                 catch
%                     message('Unable to save table and figure!');
%                 end
                
            case 'Yes'
                % Save comparison of permuted and true accuracy to CSV file
                file = fopen([fileparts(outputPath) filesep  'permAccuraciesSummary.csv'], 'w');
                
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
        clc;
        clear;
    end
end
