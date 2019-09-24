#!/bin/bash

####Notes & Comments####
howtouse() {
echo ""
echo "Register Searchlight Results to MNI"
echo "Daniel Elbich"
echo "Created: 4/15/19"
echo ""
echo ""
echo ""
echo "Usage:"
echo "sh searchlightToMNI.sh --task <text>"
echo ""
echo ""
echo " Script registers searchlight results to MNI space. If data is already"
echo " is already in MNI space, this step can be skipped."
echo ""
echo "Required arguments:"
echo ""
echo "      --path         Path to searchlight results folder"
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

#pathSearchlight='/path/to/results/from/searchlight'
prefix='./'

cd $path

subs=(./*)

for ((i=0; i<${#subs[@]}; i++))
do

sub=${subs[i]#"$prefix"}

echo $sub

cd $sub
echo "Fitting "$sub" diffusion data to MNI space..."
flirt -in *.nii -ref /opt/aci/sw/fsl/5.0.10/fsl/data/standard/MNI152_T1_2mm_brain -out "MNI_"$sub".nii.gz" -omat "MNI_"$sub".mat" -bins 256 -cost corratio -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 12  -interp trilinear

cd ..

done
