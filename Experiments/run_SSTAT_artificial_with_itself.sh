#!/bin/bash
array=$1 #3982 calculations


cd /projappl/project_2006203/MOSTA-SSTAT/Experiments
# cd csc_scratch/motifsimilarity-private/experiments/

results_path=/scratch/project_2006203/MOSTA-SSTAT

#FOXO1 FOXK1 esimerkki 624389

readarray -t filenames < /projappl/project_2006203/TFBS/PWMs_final/representatives_artificial.csv
#echo ${#filenames[@]} 24744

transfac_path=/scratch/project_2006203/TFBS/artificial_motifs/artificial_motifs_transfac/

start_ind=0
#end_ind=3982
end_ind=24744


for (( index=$start_ind; index<$end_ind; index++ )); do

i=$index
j=$index

#echo "index: "$index
echo "i: "$i
echo "j: "$j

#From i j back to index
#index=j+i*(i-1)/2


TF1=${filenames[$i]}
TF2=${filenames[$j]}

PWM1=$(basename "$TF1")
TF1=$transfac_path$PWM1
PWM1="${PWM1%.*}"

PWM2=$(basename "$TF2")
TF2=$transfac_path$PWM2
PWM2="${PWM2%.*}"

echo $TF1 > $LOCAL_SCRATCH"/"$PWM1"_"$PWM2".list"
echo $TF2 >> $LOCAL_SCRATCH"/"$PWM1"_"$PWM2".list"


if [[ $index -eq $start_ind ]]
then
     #echo $index is equal to start
     #echo $PWM1"-"$PWM2 > "../results/result_"$array".out"

     ../sstat .5 list:$LOCAL_SCRATCH"/"$PWM1"_"$PWM2".list" typeI 0.01 > $LOCAL_SCRATCH"/result_with_itself.out"



else
     echo $index is greater then start
     #echo $PWM1"-"$PWM2 >> "../results/result_"$array".out"
     ../sstat .5 list:$LOCAL_SCRATCH"/"$PWM1"_"$PWM2".list" typeI 0.01 >> $LOCAL_SCRATCH"/result_with_itself.out"


fi

rm $LOCAL_SCRATCH"/"$PWM1"_"$PWM2".list"


done

cd $LOCAL_SCRATCH
cp "result_with_itself.out" $results_path"/artificial_results_final/"

