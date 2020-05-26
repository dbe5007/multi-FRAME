%% Mask Registration
%   Editor:    Daniel Elbich
%   Created:   5/24/20
%
%   Register brain region mask to dataset for use in multivariate analyses.
%   All masks are fit to each subject separately and saved in the
%   multivaraiate parent folder
%
%   If you plan on using Freesurfer parcellations and do not already have
%   them separated, the shell script in modeling folder
%   'exportFSRegions.sh' will do so. Running that will take some time so
%   it is advisable to submit that to a batch processing system if
%   available.
%

%% Set Pipeline Parameters

% Set Analysis Parameters & Paths
% Load subject IDs, ROIs, and Condition flags
if exist('flag','var') == 0
    
    %Select parameter file is flag does not exist
    %uiopen('*.mat')
    
    [file,path]=uigetfile('*.mat','Select params file');
    filename=fullfile(path,file);
    load(filename);
    
end

%% List of Freesurfer Regions

try
    fsRegions = {'allRegions', 'fullBrainMask', 'lh_Accumbens', 'lh_Amygdala',...
        'lh_Caudate', 'lh_Hippocampus', 'lh_Insula', 'lh_Pallidum',...
        'lh_Putamen', 'lh_Substancia', 'lh_Thalamus', 'lh_WM', 'lh_bankssts',...
        'lh_caudalanteriorcingulate', 'lh_caudalmiddlefrontal',...
        'lh_corpuscallosum', 'lh_cuneus', 'lh_entorhinal', 'lh_frontalpole',...
        'lh_fusiform', 'lh_inferiorparietal', 'lh_inferiortemporal',...
        'lh_insula', 'lh_isthmuscingulate', 'lh_lateraloccipital',...
        'lh_lateralorbitofrontal', 'lh_lingual', 'lh_medialorbitofrontal',...
        'lh_middletemporal', 'lh_paracentral', 'lh_parahippocampal',...
        'lh_parsopercularis', 'lh_parsorbitalis', 'lh_parstriangularis',...
        'lh_pericalcarine', 'lh_postcentral', 'lh_posteriorcingulate',...
        'lh_precentral', 'lh_precuneus', 'lh_rostralanteriorcingulate',...
        'lh_rostralmiddlefrontal', 'lh_superiorfrontal', 'lh_superiorparietal',...
        'lh_superiortemporal', 'lh_supramarginal', 'lh_temporalpole',...
        'lh_transversetemporal', 'rh_Accumbens', 'rh_Amygdala', 'rh_Caudate',...
        'rh_Hippocampus', 'rh_Insula', 'rh_Pallidum', 'rh_Putamen',...
        'rh_Substancia', 'rh_Thalamus', 'rh_WM', 'rh_bankssts',...
        'rh_caudalanteriorcingulate', 'rh_caudalmiddlefrontal',...
        'rh_corpuscallosum', 'rh_cuneus', 'rh_entorhinal', 'rh_frontalpole',...
        'rh_fusiform', 'rh_inferiorparietal', 'rh_inferiortemporal',...
        'rh_insula', 'rh_isthmuscingulate', 'rh_lateraloccipital',...
        'rh_lateralorbitofrontal', 'rh_lingual', 'rh_medialorbitofrontal',...
        'rh_middletemporal', 'rh_paracentral', 'rh_parahippocampal',...
        'rh_parsopercularis', 'rh_parsorbitalis', 'rh_parstriangularis',...
        'rh_pericalcarine', 'rh_postcentral', 'rh_posteriorcingulate',...
        'rh_precentral', 'rh_precuneus', 'rh_rostralanteriorcingulate',...
        'rh_rostralmiddlefrontal', 'rh_superiorfrontal', 'rh_superiorparietal',...
        'rh_superiortemporal', 'rh_supramarginal', 'rh_temporalpole',...
        'rh_transversetemporal'};
end

%% Main Code

% Setup mask options
% Account for RT/regress out. No is default
maskDirFlag = {'A prior mask directory',...
    'Freesurfer (used package shell script)','I do not have one'};

[index,tf] = listdlg('PromptString','Select Mask Directory:',...
    'SelectionMode','Single','ListString',maskDirFlag,'ListSize',[280,100]);

% Confirm need for mask registration
resliceFlag=questdlg('Are masks already registered to functional data/space?',...
    'Confirm Reslicing','Yes','No','Cancel','No');

switch index
    
    case 1 % User made directory
        
        % Path to mask directory
        maskDir = inputdlg(['Enter path to mask folder data starting from: '...
            directory.Project],'File Path',[1 35],{'e.g. project/masks'});
        maskDir = [directory.Project filesep maskDir];
        
        % Get all masks in directory
        masks=dir([maskDir '/*.ni*']);
        
        % Drop-down list of found masks - select all that apply
        [index,tf] = listdlg('Name','Available Masks',...
            'PromptString','Ctrl+Click to choose masks:',...
            'ListString',masks.Name,...
            'ListSize',[280,300]);
        maskList=maskList(index);
        regions=maskName(index);
        
        for i=1:length(regions)
            datainput=maskList{i};
            reference=dataDir;
            output=[directory.Analysis filesep 'masks' filesep,roiPrefix,regions{i}];
            
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
        
    case 2 % Freesurfer ROIs parcellated via package shell script
        
        % Get Freesurfer subject directory
        fsSubjDir = inputdlg(['Enter path to freesurfer subjects directory starting after: '...
            directory.Project],'File Path',[1 35],{'e.g. path/to/subjects_dir'});
        fsSubjDir = [directory.Project filesep fsSubjDir{:}];
        
        fsRegions = dir([fsSubjDir '/' subjects{1} '/func/segmentationSplit/*.nii.gz']);
        for i=1:length(fsRegions)
            maskList{i}=fsRegions(i).name;
        end
        
        % Drop-down list of found masks - select all that apply
        [index,tf] = listdlg('Name','Available Masks',...
            'PromptString','Ctrl+Click to choose masks:',...
            'ListString',maskList,...
            'ListSize',[280,300]);
        regions=fsRegions(index);
        
        % Preallocate final variable
        %finalMasks=cell(1,length(regions));
        
        for i=1:length(subjects)
            for j=1:length(regions)
                datainput=[fsSubjDir '/' subjects{i} '/func/segmentationSplit/'...
                    regions(j).name];
                output=[directory.Analysis filesep 'masks' filesep ...
                    file(1:end-4) filesep subjects{i} filesep regions(j).name '.nii.gz'];
                
                % Set system/terminal variables
                setenv('outDir',[directory.Analysis filesep 'masks' filesep ...
                    file(1:end-4) filesep subjects{i}]);
                setenv('datainput',datainput);
                
                % Copy to mask directory
                !mkdir -p $outDir
                !cp $datainput $outDir
                
                %finalMasks{1,j}=regions{j};
            end
        end
        
    case 3
        
        msgbox(['No mask directory selected. Please create mask directory'...
            ' and rerun script!']);
        return
        
end

