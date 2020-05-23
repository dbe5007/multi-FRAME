%% Preprocess Functional Data
%  Daniel Elbich
%  Cognitive Aging & Neuroimaging Lab
%  5/22/19

% This script organizes data for preprocessing. The pipeline can accept
% both SPM and fMRIPrep input.
%
% Here, the script offers a decision point: if
% data is processed via fMRIPrep, raw data will be copied to a working
% directory for processing and motion parameters will be transformed to fit
% SPM convention (i.e. degrees and radians for translation and rotation,
% respectively).
%
% If SPM12 is the preferredn processing, the script will batch multiple 
% subjects through a preprocessing pipeline designed to collect ALL 
% AVAIABLE FUNCTIONAL RUNS USING WILDCARDS, and PREPROCESS THEM 
% ALLTOGETHER. Function scans are realigned to the first image of each run,
% respectively.
%
% This file structure is setup for BIDS compliant organization. It does not
% require that file names be in precise BIDS formatting, but does assume
% separate anat & func directories, and all subject names begin with
% 'sub-'. If this is not the case for your data, consider using this file
% structure or edit 'User Input Step 1' to fit your setup.
%
% In additon, this script will create & copy data to a separate folder for
% preprocessing. This is performed to protect integrity of raw data and
% prevent overwriting.

%% Preprocessing Decision Point

if exist('commandFlag','var')==0
    preprocPipe = questdlg('Select Preprocessing Pipeline',...
        'Confirm Preprocessing Pipeline',...
        'fMRIPrep','SPM12','Cancel','SPM12');
else
    preprocPipe = 'MVPA';
end

%% Initial Path & File Setup

% Project directories
projDir = '/director/to/project';

switch preprocPipe
    case 'SPM12'
        
        % Add SPM12 to path
        addpath('/directory/to/spm12');
        dataDir = [projDir '/rawdata'];
        processDir = [projDir '/derivatives/spmPreprocessing'];
        setenv('dataDir',dataDir);
        setenv('processDir',processDir);
        
        % Set file and folder wildcard expressions
        % Specify a regular expression (google regular expressions) that will
        % select only the raw image functional & anatomical image respectively.
        
        wildcard.func = '^*bold.nii';
        wildcard.anat = '^*T1w.nii';
        wildcard.run = 'run.*';
        
        % Set output directories for preprocessed data
        directories.func    = processDir;
        directories.anat    = processDir;
        directories.psfiles = [processDir '/psfiles'];
        !mkdir -p $processDir/psfiles;
        
    case 'fMRIPrep'
        dataDir = [projDir '/derivatives/fmriprep'];
end

if exist('commandFlag','var')==0
    
    %% Subjects
    % Remove non-subject directories
    subjListFlag = questdlg('Choose subjects to process:', ...
        'Subject Processing', ...
        'Automated','User selected','Cancel','Automated');
    
    % Handle response
    switch subjListFlag
        case 'Automated'
            automationFlag = questdlg('Choose automation process:', ...
                'Subject Processing', ...
                'All Subjects','Only Unprocessed','Cancel','All Subjects');
            
            switch automationFlag
                case 'Cancel'
                    clear;
                    disp('User cancelled process...');
                    return;
            end
            
            funcFolders = dir(dataDir);
            
            for i=1:length(funcFolders)
                subFlag(i) = ~isempty(strfind(funcFolders(i).name,'sub-'));
            end
            funcFolders = funcFolders(subFlag);
            
            % Remove errant files
            dirFlags = [funcFolders.isdir];
            subfolders = funcFolders(dirFlags);
            
            for i=1:length(subfolders)
                subjects{i} = subfolders(i).name;
            end
            
        case 'User selected'
            subjects = inputdlg('Enter subject ID. For multiple separate with space:'...
                ,'Subject ID Input',[1 50]);
            subjects = split(subjects,[' ',"'"]);
            
            % Whitespace check
            for i=1:length(subjects)
                remove(i) = ~isempty(subjects{i});
            end
            subjects = subjects(remove);
            
        case 'Cancel'
            clear;
            disp('User cancelled process...');
            return;
    end
    
else
    
    subjListFlag = 'All Subjects';
    funcFolders = dir(dataDir);
    
    for i=1:length(funcFolders)
        subFlag(i) = ~isempty(strfind(funcFolders(i).name,'sub-'));
    end
    funcFolders = funcFolders(subFlag);
    
    % Remove errant files
    dirFlags = [funcFolders.isdir];
    subfolders = funcFolders(dirFlags);
    
    for i=1:length(subfolders)
        subjects{i} = subfolders(i).name;
    end
    
end

%% Main Code

switch preprocPipe
    case 'SPM12'
        % Initialize SPM
        spm('defaults', 'FMRI'); % load SPM default options
        spm_jobman('initcfg')    % Configure the SPM job manger
        
        % Loop for all subjects
        for csub = subjects
            
            switch subjListFlag
                case 'Only Unprocessed'
                    if exist([processDir filesep csub{:} filesep...
                            'func/ashburnerReferenceImage.nii.gz'],'file')==2
                        continue;
                    end
            end
            
            % Copy data from raw BIDS dir
            setenv('subject',csub{1,1});
            
            % Copy & unzip T1
            !mkdir -p $processDir/$subject/anat
            !cp $dataDir/$subject/anat/sub* $processDir/$subject/anat
            !gunzip $processDir/$subject/anat/sub*
            
            % Copy & unzip functionals
            !mkdir -p $processDir/$subject/func
            !cp $dataDir/$subject/func/*nii.gz $processDir/$subject/func
            !gunzip $processDir/$subject/func/sub*
            
            % Create separate directories for each run
            !mkdir -p $processDir/$subject/func/{run1,run2,run3,run4}
            !mv $processDir/$subject/func/*run-1* $processDir/$subject/func/run1
            !mv $processDir/$subject/func/*run-2* $processDir/$subject/func/run2
            !mv $processDir/$subject/func/*run-3* $processDir/$subject/func/run3
            !mv $processDir/$subject/func/*run-4* $processDir/$subject/func/run4
            
            % Create the path to this subjects' functional folder
            subjFunc = fullfile(directories.func, csub{:}, 'func');
            
            % Select run folders
            runs = cellstr(spm_select('FPList', subjFunc, 'dir', wildcard.run));
            
            %% Defining Session Parameters
            % Get paths for each volume of 4D nifti images for each run
            for curRun = 1:length(runs)
                % Pull all images/volumes into single variable
                runDir  = runs{curRun};
                images.raw{curRun} = cellstr(spm_select('ExtFPList', runDir, wildcard.func, Inf));
            end
            
            % Set number of preprocessing parameters
            params = 7;
            matlabbatch=cell(1,7);
            
            for i=1:params
                switch i
                    case 1
                        % Run Independent Realightment Parameters
                        matlabbatch{i}.spm.spatial.realign.estwrite.data = images.raw; % Paths to all images to realign for this run
                        matlabbatch{i}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
                        matlabbatch{i}.spm.spatial.realign.estwrite.eoptions.sep     = 4;
                        matlabbatch{i}.spm.spatial.realign.estwrite.eoptions.fwhm    = 5;
                        matlabbatch{i}.spm.spatial.realign.estwrite.eoptions.rtm     = 1; % Register to mean of the images
                        %matlabbatch{i}.spm.spatial.realign.estwrite.eoptions.rtm     = 0; % Register to 1st image in series
                        matlabbatch{i}.spm.spatial.realign.estwrite.eoptions.interp  = 2;
                        matlabbatch{i}.spm.spatial.realign.estwrite.eoptions.wrap    = [0 0 0];
                        matlabbatch{i}.spm.spatial.realign.estwrite.eoptions.weight  = '';
                        %matlabbatch{i}.spm.spatial.realign.estwrite.roptions.which   = [0 1]; % Create Mean image only
                        matlabbatch{i}.spm.spatial.realign.estwrite.roptions.which   = [2 1]; % Create resliced images of all runs and mean
                        matlabbatch{i}.spm.spatial.realign.estwrite.roptions.interp  = 4;
                        matlabbatch{i}.spm.spatial.realign.estwrite.roptions.wrap    = [0 0 0];
                        matlabbatch{i}.spm.spatial.realign.estwrite.roptions.mask    = 1;
                        matlabbatch{i}.spm.spatial.realign.estwrite.roptions.prefix  = 'r';
                        
                    case 2
                        % Run Independent Slicetiming Parameters
                        for a=1:length(runs)
                            images.reslice{1,a} = strrep(images.raw{1,a},[csub{:} '_task'],...
                                ['r' csub{:} '_task']);
                        end
                        
                        matlabbatch{i}.spm.temporal.st.scans = images.reslice;
                        matlabbatch{i}.spm.temporal.st.nslices  = 58;
                        matlabbatch{i}.spm.temporal.st.tr       = 2.5;
                        matlabbatch{i}.spm.temporal.st.ta       = 2.5-(2.5/58);
                        %matlabbatch{i}.spm.temporal.st.so       = [2:2:58 1:2:58];   %Even slices collected 1st, then odds; start foot to head
                        matlabbatch{i}.spm.temporal.st.so       = [1:2:58 2:2:58];   %Odd slices collected 1st, then odds; start foot to head
                        matlabbatch{i}.spm.temporal.st.refslice = 2;
                        matlabbatch{i}.spm.temporal.st.prefix   = 'a';
                        
                    case 3
                        % Ashburner Fix for Better Normalization
                        meanImage = dir([subjFunc '/run1/sub*']);
                        matlabbatch{i}.spm.util.imcalc.input(1)       = {[meanImage(1).folder '/mean' meanImage(1).name]};
                        matlabbatch{i}.spm.util.imcalc.output         = 'ashburnerReferenceImage';
                        matlabbatch{i}.spm.util.imcalc.outdir         = {subjFunc};
                        matlabbatch{i}.spm.util.imcalc.expression     = 'i1 + randn(size(i1))*50';
                        matlabbatch{i}.spm.util.imcalc.var            = struct('name', {}, 'value', {});
                        matlabbatch{i}.spm.util.imcalc.options.dmtx   = 0;
                        matlabbatch{i}.spm.util.imcalc.options.mask   = 0;
                        matlabbatch{i}.spm.util.imcalc.options.interp = 1;
                        matlabbatch{i}.spm.util.imcalc.options.dtype  = 4;
                        
                    case 4
                        % Run Independent Coregistration Parameters
                        % Ashburner
                        anatDir = fullfile(directories.anat, csub{:}, 'anat');
                        matlabbatch{i}.spm.spatial.coreg.estimate.ref = {spm_select('ExtFPListRec', anatDir, wildcard.anat)};
                        matlabbatch{i}.spm.spatial.coreg.estimate.source = {[processDir ...
                            filesep csub{:} filesep 'func' filesep ...
                            matlabbatch{i-1}.spm.util.imcalc.output '.nii']};
                        
                        for a=1:length(runs)
                            images.slicetime{a,1} = [runs{a} '/ar' csub{:} '_task-mp_run-' num2str(a) '_bold.nii'];
                        end
                        
                        matlabbatch{i}.spm.spatial.coreg.estimate.other = images.slicetime;
                        matlabbatch{i}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
                        matlabbatch{i}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
                        matlabbatch{i}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                        matlabbatch{i}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
                        
                    case 5
                        % Run Independent Segmentation Parameters
                        spm_segment_image = which('TPM.nii'); % find the TPM image on the MATLAB search path
                        matlabbatch{i}.spm.spatial.preproc.channel.vols(1)  = matlabbatch{i-1}.spm.spatial.coreg.estimate.ref;
                        matlabbatch{i}.spm.spatial.preproc.channel.biasreg  = 0.001;
                        matlabbatch{i}.spm.spatial.preproc.channel.biasfwhm = 60;
                        matlabbatch{i}.spm.spatial.preproc.channel.write    = [0 0];
                        matlabbatch{i}.spm.spatial.preproc.tissue(1).tpm    = {[spm_segment_image ',1']};
                        matlabbatch{i}.spm.spatial.preproc.tissue(1).ngaus  = 1;
                        matlabbatch{i}.spm.spatial.preproc.tissue(1).native = [1 0];
                        matlabbatch{i}.spm.spatial.preproc.tissue(1).warped = [0 0];
                        matlabbatch{i}.spm.spatial.preproc.tissue(2).tpm    = {[spm_segment_image ',2']};
                        matlabbatch{i}.spm.spatial.preproc.tissue(2).ngaus  = 1;
                        matlabbatch{i}.spm.spatial.preproc.tissue(2).native = [1 0];
                        matlabbatch{i}.spm.spatial.preproc.tissue(2).warped = [0 0];
                        matlabbatch{i}.spm.spatial.preproc.tissue(3).tpm    = {[spm_segment_image ',3']};
                        matlabbatch{i}.spm.spatial.preproc.tissue(3).ngaus  = 2;
                        matlabbatch{i}.spm.spatial.preproc.tissue(3).native = [1 0];
                        matlabbatch{i}.spm.spatial.preproc.tissue(3).warped = [0 0];
                        matlabbatch{i}.spm.spatial.preproc.tissue(4).tpm    = {[spm_segment_image ',4']};
                        matlabbatch{i}.spm.spatial.preproc.tissue(4).ngaus  = 3;
                        matlabbatch{i}.spm.spatial.preproc.tissue(4).native = [1 0];
                        matlabbatch{i}.spm.spatial.preproc.tissue(4).warped = [0 0];
                        matlabbatch{i}.spm.spatial.preproc.tissue(5).tpm    = {[spm_segment_image ',5']};
                        matlabbatch{i}.spm.spatial.preproc.tissue(5).ngaus  = 4;
                        matlabbatch{i}.spm.spatial.preproc.tissue(5).native = [1 0];
                        matlabbatch{i}.spm.spatial.preproc.tissue(5).warped = [0 0];
                        matlabbatch{i}.spm.spatial.preproc.tissue(6).tpm    = {[spm_segment_image ',6']};
                        matlabbatch{i}.spm.spatial.preproc.tissue(6).ngaus  = 2;
                        matlabbatch{i}.spm.spatial.preproc.tissue(6).native = [0 0];
                        matlabbatch{i}.spm.spatial.preproc.tissue(6).warped = [0 0];
                        matlabbatch{i}.spm.spatial.preproc.warp.mrf         = 1;
                        matlabbatch{i}.spm.spatial.preproc.warp.cleanup     = 1;
                        matlabbatch{i}.spm.spatial.preproc.warp.reg         = [0 0.001 0.5 0.05 0.2];
                        matlabbatch{i}.spm.spatial.preproc.warp.affreg      = 'mni';
                        matlabbatch{i}.spm.spatial.preproc.warp.fwhm        = 0;
                        matlabbatch{i}.spm.spatial.preproc.warp.samp        = 3;
                        matlabbatch{i}.spm.spatial.preproc.warp.write       = [0 1];
                        
                    case 6
                        % Run Indepenedent Normalization: Normalizaing the Functional Images
                        count = 1;
                        for a=1:length(runs)
                            for j=1:length(images.raw{a})
                                images.norm{count,1}=images.reslice{a}{j};
                                images.norm{count,1}=strrep(images.norm{count,1},...
                                    ['r' csub{:} '_task'],['ar' csub{:} '_task']);
                                count=count+1;
                            end
                        end
                        
                        images.norm{count,1} = matlabbatch{i-2}.spm.spatial.coreg.estimate.ref{1};
                        matlabbatch{i}.spm.spatial.normalise.estwrite.subj.vol = {[anatDir '/' csub{:} '_T1w.nii']};
                        matlabbatch{i}.spm.spatial.normalise.estwrite.subj.resample = images.norm;
                        matlabbatch{i}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
                        matlabbatch{i}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
                        matlabbatch{i}.spm.spatial.normalise.estwrite.eoptions.tpm = {'path/to/spm12/tpm/TPM.nii'};
                        matlabbatch{i}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
                        matlabbatch{i}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
                        matlabbatch{i}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
                        matlabbatch{i}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
                        matlabbatch{i}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
                            78 76 85];
                        matlabbatch{i}.spm.spatial.normalise.estwrite.woptions.vox = [3 3 2];
                        matlabbatch{i}.spm.spatial.normalise.estwrite.woptions.interp = 4;
                        matlabbatch{i}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';
                        
                        
                    case 7
                        %Run Independent Smoothing Parameters
                        for a=1:(count-1)
                            images.smooth{a,1}=strrep(images.norm{a,1},...
                                ['ar' csub{:} '_task'],['war' csub{:} '_task']);
                        end
                        
                        matlabbatch{i}.spm.spatial.smooth.data = images.smooth;
                        matlabbatch{i}.spm.spatial.smooth.fwhm    = [6 6 4];
                        matlabbatch{i}.spm.spatial.smooth.dtype   = 0;
                        matlabbatch{i}.spm.spatial.smooth.im      = 0;
                        matlabbatch{i}.spm.spatial.smooth.prefix  = 's';
                end
            end
            
            %% Run batch
            % Set the flag to 1 to look at the parameters interactively in the GUI (e.g. debug)
            flag = 2;
            % Configure spm graphics window. Ensures a .ps file is saved during preprocessing
            spm_figure('GetWin','Graphics');
            
            % Make psfiles the working directory. Ensures .ps file is saved in this directory
            cd(directories.psfiles);
            
            % Run preprocessing
            spm_jobman('run', matlabbatch);
            
            % Rename the ps file from "spm_CurrentDate.ps" to "SubjectID.ps"
            try
                temp = date;
                newDate = [temp(end-3:end) temp(4:6) temp(1:2)];
                movefile(['spm_' newDate '.ps'],sprintf('%s.ps',csub{:}));
            catch
                % Catch for if processing ends on different day
                % (e.g. changes in date from start to end)
                temp = date;
                newDay = str2double(temp(1:2))-1;
                if newDay <=9
                    newDate = [temp(end-3:end) temp(4:6) '0' num2str(newDay)];
                else
                    newDate = [temp(end-3:end) temp(4:6) num2str(newDay)];
                end
                movefile(['spm_' newDate '.ps'],sprintf('%s.ps',csub{:}));
                continue;
            end
            
            clear images matlabbatch newBatch;
            
            % Gzip files to save space
            !find $processDir/$subject -type f -name '*.nii' -exec gzip "{}" \;
            
        end
        
    case 'fMRIPrep'
        
        dataDir = [projDir '/fmriprep'];
        for csub = subjects
            
        end
end

