#!/bin/sh

####Notes & Comments####
help() {
echo ""
echo "Create multi-FRAME PBS Scripts"
echo "Daniel Elbich"
echo "Cogntive, Aging, and Neurogimaging Lab"
echo "Created: 5/29/19"
echo ""
echo ""
echo " Creates project specific PBS job files for command line use of multi-FRAME pipeline. Files "
echo " are saved in the project directory."
echo ""
echo ""
echo "Usage:"
echo "sh createPBSScript.sh --paramsFile <params.mat filename> --paramsDir <params.mat directory> --projectDir <Project Directory>"
echo ""
echo " Required arguments:"
echo ""
echo "	    --paramsFile    Parameter filename as output by createparams.m script"
echo "	    --paramsDir     Parameter file directory created by createparams.m script"
echo ""
echo ""
exit 1
}
[ "$1" = "--help" ] && help


#Argument check
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
--paramsFile) paramsFile="$2"
shift # past argument
shift # past value
;;
--paramsDir) paramsDir="$2"
shift # past argument
shift # past value
;;
--projectDir) projectDir="$2"
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

#####OTHER PBS FLAGS#####
#Delayed qsub submission line - military time
#PBS -a 2200

### PBS Variables ###
ALLOCATION=open
NODES=1
PPN=4
WALLTIME='48:00:00'
MEM=4gb

## Save pbs text file to directory##
FILE=$projectDir/PBS/'PBS_multi-FRAME_'$paramsFile'.txt'
OUTPUT=$projectDir/PBS/output
ERROR=$projectDir/PBS/error
JOBNAME='multi-FRAMEProc' #MAX 15 CHARACTERS!!!

#paramsFilename='paramsFilename='$paramsFile
paramsFilename=$(echo 'paramsFilename='$paramsFile'.mat')
paramsFileFullpath=$(echo 'paramsFileFullpath='$paramsDir/$paramsFile'.mat')

preproc=$(echo 'matlab -nodisplay -nosplash -nodesktop -r "commandFlag=1;file='"'"$paramsFile'.mat'"'"';load('"'"$paramsDir/$paramsFile'.mat'"'"');run(which('"'"multi-FRAME/modeling/preprocessData.m"'"'));quit"')
specify=$(echo 'matlab -nodisplay -nosplash -nodesktop -r "commandFlag=1;file='"'"$paramsFile'.mat'"'"';load('"'"$paramsDir/$paramsFile'.mat'"'"');run(which('"'"multi-FRAME/modeling/speciyModel.m"'"'));quit"')
estimate=$(echo 'matlab -nodisplay -nosplash -nodesktop -r "commandFlag=1;file='"'"$paramsFile'.mat'"'"';load('"'"$paramsDir/$paramsFile'.mat'"'"');run(which('"'"multi-FRAME/modeling/estimateModel.m"'"'));quit"')
run=$(echo 'matlab -nodisplay -nosplash -nodesktop -r "commandFlag=1;file='"'"$paramsFile'.mat'"'"';load('"'"$paramsDir/$paramsFile'.mat'"'"');run(which('"'"multi-FRAME/analyses/runMVPAClassification.m"'"'));quit"')

fullPreproc=$(echo 'matlab -nodisplay -nosplash -nodesktop -r "commandFlag=1;file='"'"$paramsFile'.mat'"'"';load('"'"$paramsDir/$paramsFile'.mat'"'"');run(which('"'"multi-FRAME/modeling/preprocessData.m"'"'));run(which('"'"multi-FRAME/modeling/speciyModel.m"'"'));run(which('"'"multi-FRAME/modeling/estimateModel.m"'"'));quit"')

/bin/cat <<EOM >$FILE
#PBS -A $ALLOCATION
#PBS -l nodes=$NODES:ppn=$PPN
#PBS -l walltime=$WALLTIME
#PBS -l pmem=$MEM
#PBS -o $OUTPUT
#PBS -e $ERROR
#PBS -N $JOBNAME

## Load MATLAB ##
module load matlab/R2017b

## Set Parameter File Information ##
$paramsFilename
$paramsFileFullpath

## Run Preprocessing ##
$preproc

## Specify Model ##
$specify

## Estimate Model ##
$estimate

## Run MVPA Classification ##
$run

EOM


