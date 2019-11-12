#!/bin/bash

####Notes & Comments####
howtouse() {
echo ""
echo "One Sample T-test for Searchlight Results"
echo "Daniel Elbich"
echo "Created: 4/15/19"
echo ""
echo ""
echo ""
echo "Usage:"
echo "sh facesceneSearchlightOneSample.sh --task <text>"
echo ""
echo ""
echo " Script concatenates all subject searchlight results into 4D nifti file. Results"
echo " are then entered into one sample t-test using FSL randomise. Follows procedure "
echo " outlined in Bowman, Chamberlain, & Dennis (2019)."
echo ""
echo "Required arguments:"
echo ""
echo "      --path         Flag for task (enc = Encoding, ret = Retrieval)"
echo ""
echo ""
exit 1
}
[ "$1" = "--help" ] && howtouse

#Argument check
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
--task) path="$2"
shift # past argument
shift # past value
;;
*)    # unknown option
POSITIONAL+=("$1") # save it in an array for later
shift # past argument
;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

echo $path

#Change directory to analysis directory#
cd $path

#Get list of all subjects#
subs=(*)

#Copy data to folder to concatenate all
mkdir Concat

##Loop through all subject directories and copy data to Concat folder
for ((i=0; i<${#subs[@]}; i++))
do

#sub=${subs[i]#"$prefix"}
sub=${subs[i]}

echo $sub

#cp $sub/*.nii Concat

##Create binarized mask from results##
fslmaths $sub/*.nii -bin $sub/$sub"_bin_mask.nii.gz"

##Threshold map to 0.5 for randomise - sets chance to 50% - and threshold with binarized results##
fslmaths $sub/*.nii -sub 0.5 $sub/$sub"_thresholded_0.5.nii.gz"

#Mask with with binarized results#
fslmaths $sub/$sub"_thresholded_0.5.nii.gz" $sub/$sub"_bin_mask.nii.gz" $sub/$sub"_thresholded_0.5_masked.nii.gz"

#Fix negative 0 issue (i.e. -0) in non-searchlight voxels#
fslmaths $sub/$sub"_thresholded_0.5_masked.nii.gz" -add 0 $sub/$sub"_thresholded_0.5_masked.nii.gz"

#Copy to Concat Folder
cp $sub/$sub"_thresholded_0.5_masked.nii.gz" Concat/$sub"_thresholded_0.5.nii"

done

#Make output directory in Concat folder#
mkdir $path/Concat/Automatic
fslmerge -a $path/Concat/Automatic/"searchlightConcatenated.nii.gz" $path/Concat/*.nii.gz
randomise -i $path/Concat/Automatic/"searchlightConcatenated.nii.gz" -o $path/Concat/Automatic/"oneSampT" -1 -v 2 -T

