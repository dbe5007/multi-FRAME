#!/bin/bash

####Notes & Comments####
help() {
echo ""
echo "Create Freesurfer Masks"
echo "Daniel Elbich"
echo "Cogntive, Aging, and Neurogimaging Lab"
echo "Created: 4/2/19"
echo ""
echo ""
echo " Splits the Freesurfer cortical and subcortical segmentation into separate mask files. The"
echo " input is the transformed & registered 'aparcaseg' file in the 'func' subject directory."
echo " All masks are output as gzip NiFTi files for use in other packages (e.g. FSL, SPM). Requires "
echo " FSL be in the path, as well as the BIDS data output hierarchy created via from fMRIPrep."
echo ""
echo ""
echo "Usage:"
echo "sh splitFMRIPrepFSMasks.sh --subjList <subjectIDs> --task <task>"
echo ""
echo " Required arguments:"
echo ""
echo "      --subjList      Single subject ID or text file containing list of subject IDs"
echo "      --task          Task name for correct aparc+aseg file selection. Still required if only 1 task"
echo ""
echo " Optional arguments (You may optionally specify one or more of): "
echo ""
echo "	    --fsSubjDir     Freesurfer subjects directory (change if you do not"
echo "                          to use the default setup from sourcing Freesurfer)"
echo "      --multipleruns  Include if there are multiple runs of a single task. Will default to using 1st run"
echo ""
echo ""
exit 1
}
[ "$1" = "--help" ] && help

#Arguement check
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in

--subjList) subjList="$2"
shift # past argument
shift # past value
;;
--task) task="$2"
shift # past argument
shift # past value
;;
--fsSubjDir) fsSubjDir="$2"
shift # past argument
shift # past value
;;
--multipleruns) multipleruns="$1"
shift # past argument
;;
*)    # unknown option
POSITIONAL+=("$1") # save it in an array for later
shift # past argument
;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

###Read Subject IDs###
if [[ -f $subjList ]]; then
	i=0
	while read -r LINE || [[ -n $LINE ]]; do
		subs[i]=$LINE
		let "i++"
	done < $subjList
else
	subs=$subjList
fi

##Load in freesurfer and fsl - FOR SERVER USE ONLY!!!##
module load fsl
module load freesurfer
source /opt/aci/sw/freesurfer/6.0.0/SetUpFreeSurfer.sh

#Setup freesurfer subject directory##
if [ -z ${fsSubjDir+x} ]
then
	echo "Using default subjects directory..."
else
	export SUBJECTS_DIR=$fsSubjDir
	echo $SUBJECTS_DIR
fi
echo "Subjects directory is: " $SUBJECTS_DIR

##Color LUTs for Aseg parcellation - refer to freesurfer webpage##
#https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/AnatomicalROI/FreeSurferColorLUT#
# Includes subcortical regions excluding all ventricles #
#subcortColorCode=(9 11 12 13 17 18 19 26 27 48 50 51 52 53 54 55 58 59)
#subcortColorName=('lh_Thalamus' 'lh_Caudate' 'lh_Putamen' 'lh_Pallidum' 'lh_Hippocampus' 'lh_Amygdala' 'lh_Insula' 'lh_Accumbens' 'lh_Substancia' 'rh_Thalamus' 'rh_Caudate' 'rh_Putamen' 'rh_Pallidum' 'rh_Hippocampus' 'rh_Amygdala' 'rh_Insula' 'rh_Accumbens' 'rh_Substancia')

subcortColorCode=(6 7 8 9 10 11 12 13 17 18 19 20 26 27 28 45 46 47 48 49 50 51 52 53 54 55 56 58 59 60)
subcortColorName=('lh_Cerebellum_Ext' 'lh_Cerebellum_WM' 'lh_Cerebellum_Cortex' 'lh_Thalamus' 'lh_Thalamus_Proper' 'lh_Caudate' 'lh_Putamen' 'lh_Pallidum' 'lh_Hippocampus' 'lh_Amygdala' 'lh_Insula' 'lh_Operculum' 'lh_Accumbens' 'lh_Substancia' 'lh_VentralDC' 'rh_Cerebellum_Ext' 'rh_Cerebellum_WM' 'rh_Cerebellum_Cortex' 'rh_Thalamus' 'rh_Thalamus_Proper' 'rh_Caudate' 'rh_Putamen' 'rh_Pallidum' 'rh_Hippocampus' 'rh_Amygdala' 'rh_Insula' 'rh_Operculum' 'rh_Accumbens' 'rh_Substancia' 'rh_VentralDC')

# Cortical labels #
aparcColorName=('unknown' 'bankssts' 'caudalanteriorcingulate' 'caudalmiddlefrontal' 'corpuscallosum' 'cuneus' 'entorhinal' 'fusiform' 'inferiorparietal' 'inferiortemporal' 'isthmuscingulate' 'lateraloccipital' 'lateralorbitofrontal' 'lingual' 'medialorbitofrontal' 'middletemporal' 'parahippocampal' 'paracentral' 'parsopercularis' 'parsorbitalis' 'parstriangularis' 'pericalcarine' 'postcentral' 'posteriorcingulate' 'precentral' 'precuneus' 'rostralanteriorcingulate' 'rostralmiddlefrontal' 'superiorfrontal' 'superiorparietal' 'superiortemporal' 'supramarginal' 'frontalpole' 'temporalpole' 'transversetemporal' 'insula')

# Multiple run flag #
if [ "$multipleruns" = "--multipleruns" ]
then
	task=$task"_run-1"
fi

## Change to project directory and get subject list ##
cd $SUBJECTS_DIR

## Loop for all subjects ##
for ((i=0; i<${#subs[@]}; i++))
do

## Change to project directory and get subject list ##
#cd $SUBJECTS_DIR

sub=${subs[i]}
echo $sub

## Get list of labels and make output directory for conversion ##
mkdir -p $SUBJECTS_DIR/$sub/func/segmentationSplit #fMRIPrep Pre-Registered APARC+ASEG Directory

## Convert aparc regions to nii ##
lhCount=1000  #Left hemisphere regions range from 1000 to 1035
rhCount=2000  #Right hemisphere regions range from 2000 to 2035

### WM Masks ###
## Thresholds map to specific color value and saves as nifti ###
fslmaths $SUBJECTS_DIR/$sub/func/*$task*aparc*.nii.gz -uthr 41 -thr 41 $SUBJECTS_DIR/$sub/func/segmentationSplit/'lh_WM.nii.gz'
fslmaths $SUBJECTS_DIR/$sub/func/*$task*aparc*.nii.gz -uthr 2 -thr 2 $SUBJECTS_DIR/$sub/func/segmentationSplit/'rh_WM.nii.gz'

## Binarizes voxels using color value ##
fslmaths $SUBJECTS_DIR/$sub/func/segmentationSplit/'lh_WM.nii.gz' -div 41 $SUBJECTS_DIR/$sub/func/segmentationSplit/'lh_WM.nii.gz'
fslmaths $SUBJECTS_DIR/$sub/func/segmentationSplit/'rh_WM.nii.gz' -div 2 $SUBJECTS_DIR/$sub/func/segmentationSplit/'rh_WM.nii.gz'

## Unzip files to separate from inclusion in grey matter mask creation (see below) ##
gunzip $SUBJECTS_DIR/$sub/func/segmentationSplit/lh_WM.nii.gz
gunzip $SUBJECTS_DIR/$sub/func/segmentationSplit/rh_WM.nii.gz


echo "Converting cortical parcellations..."
for ((ii=0; ii<${#aparcColorName[@]}; ii++))
do

## Thresholds map to specific color value and saves as nifti ###
fslmaths $SUBJECTS_DIR/$sub/func/*$task*aparc*.nii.gz -uthr $lhCount -thr $lhCount $SUBJECTS_DIR/$sub/func/segmentationSplit/'lh_'${aparcColorName[ii]}'.nii.gz'
fslmaths $SUBJECTS_DIR/$sub/func/*$task*aparc*.nii.gz -uthr $rhCount -thr $rhCount $SUBJECTS_DIR/$sub/func/segmentationSplit/'rh_'${aparcColorName[ii]}'.nii.gz'

## Binarizes voxels using color value ##
fslmaths $SUBJECTS_DIR/$sub/func/segmentationSplit/'lh_'${aparcColorName[ii]}'.nii.gz' -div $lhCount $SUBJECTS_DIR/$sub/func/segmentationSplit/'lh_'${aparcColorName[ii]}'.nii.gz'
fslmaths $SUBJECTS_DIR/$sub/func/segmentationSplit/'rh_'${aparcColorName[ii]}'.nii.gz' -div $rhCount $SUBJECTS_DIR/$sub/func/segmentationSplit/'rh_'${aparcColorName[ii]}'.nii.gz'

let "lhCount++"
let "rhCount++"

done


## Convert subcortical regions to nii ##
echo "Converting subcortical parcellations..."
for ((ii=0; ii<${#subcortColorCode[@]}; ii++))
do

region=${subcortColorName[ii]}
#echo $region

# Thresholds map to specific color value and saves as nifti #
fslmaths $SUBJECTS_DIR/$sub/func/*$task*aparc*.nii.gz -uthr ${subcortColorCode[ii]} -thr ${subcortColorCode[ii]} $SUBJECTS_DIR/$sub/func/segmentationSplit/$region'.nii.gz'

# Binarizes voxels using color value #
fslmaths $SUBJECTS_DIR/$sub/func/segmentationSplit/$region'.nii.gz' -div ${subcortColorCode[ii]} $SUBJECTS_DIR/$sub/func/segmentationSplit/$region'.nii.gz'

done

## Create graymatter mask ##
cd $SUBJECTS_DIR/$sub/func/segmentationSplit
#echo $PWD

# List for all mask files #
masks=(*.nii.gz)

# Iteratively adds all mask files into single file #
echo "Creating greymatter mask..."
for ((j=0; j<${#masks[@]}; j++))
do

if [ $j == 0 ]; then
	cp ${masks[$j]} allRegions.nii.gz
else
	fslmaths allRegions.nii.gz -add ${masks[$j]} allRegions.nii.gz

fi

done

echo "Create full brain mask..."
fslmaths allRegions.nii.gz -add lh_WM.nii -add rh_WM.nii fullBrainMask.nii.gz
fslmaths fullBrainMask.nii.gz -bin fullBrainMask.nii.gz

gzip rh_WM.nii
gzip lh_WM.nii

done


