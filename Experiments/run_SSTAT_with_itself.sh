#!/bin/bash
array=$1 #3982 calculations


cd /projappl/project_2006203/MOSTA-SSTAT/Experiments
# cd csc_scratch/motifsimilarity-private/experiments/

results_path=/scratch/project_2006203/MOSTA-SSTAT

#FOXO1 FOXK1 esimerkki 624389

start_ind=0
#end_ind=3982
end_ind=3294


for (( index=$start_ind; index<$end_ind; index++ )); do

i=$index
j=$index

#echo "index: "$index
echo "i: "$i
echo "j: "$j

#From i j back to index
#index=j+i*(i-1)/2


#filenames=( $(cut -d ',' -f1 ../../TFBS/Results/tomtom/tomtom/filenames.csv ) )
filenames=( $(cut -d ',' -f1 ../../TFBS/PWMs_final/filenames.csv ) )
length=${#filenames[@]}
TF_PATH=../../TFBS/
TF1=${filenames[$i]}
TF2=${filenames[$j]}

#convert "pwms" to "transfac"
TF1="${TF1/pwms/transfac}"
TF2="${TF2/pwms/transfac}"

#Yimengs motifs
#TF1="${TF1/pfm_from_newData/transfac}"
#TF2="${TF2/pfm_from_newData/transfac}"


#The .pfm files  need to be tab separated, if it is already tab separated do nothing
#Now they are all tab separated so no need for this
#awk -v OFS="\t" '$1=$1' $TF_PATH$TF1 > $TF_PATH${TF1%.*}"_tab.pfm"
#awk -v OFS="\t" '$1=$1' $TF_PATH$TF2 > $TF_PATH${TF2%.*}"_tab.pfm"


echo $TF1
echo $TF2
#head -n 1 $TF_PATH$TF1 | awk '{print NF}'
#10
#head -n 1 $TF_PATH$TF2 | awk '{print NF}'
motifnames=( $(cut -d ',' -f1 ../../TFBS/PWMs_final/motifnames.csv ) )
#motifnames=( $(cut -d ',' -f1 ../../TFBS/Results/tomtom/tomtom/motifnames.csv ) )

PWM1=${motifnames[$i]}
PWM2=${motifnames[$j]}

#PWM1=${TF1#*${TF1%/*}}
#PWM1=${PWM1#*/}
#PWM1=${PWM1%.*}

#PWM2=${TF2#*${TF2%/*}}
#PWM2=${PWM2#*/}
#PWM2=${PWM2%.*}

#we are in motifsimilarity/experiments/slurm

#If you use > instead of >>, you will overwrite destfile rather than add to it.

echo "../../TFBS/"$TF1 > $results_path/data/$PWM1"_"$PWM2".list"
echo "../../TFBS/"$TF2 >> $results_path/data/$PWM1"_"$PWM2".list"

if [[ $index -eq $start_ind ]]
then
     echo $index is equal to start
     #echo $PWM1"-"$PWM2 > "../results/result_"$array".out"

     ../sstat .5 list:$results_path/data/$PWM1"_"$PWM2".list" typeI 0.01 > $results_path"/results_final/result_with_itself.out"



else
     echo $index is greater then start
     #echo $PWM1"-"$PWM2 >> "../results/result_"$array".out"
     ../sstat .5 list:$results_path/data/$PWM1"_"$PWM2".list" typeI 0.01 >> $results_path"/results_final/result_with_itself.out"


fi

rm $results_path/data/$PWM1"_"$PWM2".list"


done
#rm -rf errs && mkdir -p errs
#rm -rf outs && mkdir -p outs

#Also remove earlier results if needed
#rm -rf <results_path> && mkdir -p <results_path>

#n=3982

#for(index in seq(0,5,1) ){

  #i=floor( (1+sqrt(1+8*index) )/2)
  #j=index-i*(i-1)/2
  #print(index)
  #print(1+sqrt(1+8*index))
  #print(paste0("i: ", i))
  #print(paste0("j: ", j))

#}

#for (( index=0; index<6; index++ )); do

	#tmp1=$((1+8*$index))
	#tmp=`echo sqrt "($tmp1)" | bc -l`
	#tmp1=`echo "($tmp)+1" | bc -l`
	#i=`echo "($tmp1)/2" | bc -l`
	#i=${i%.*}
 	#j=$(($index-$i*($i-1)/2))

	#a=$(awk -v x=$tmp 'BEGIN{print 1+sqrt(x)}')
	#echo "index: "$index
	#echo "i: "$i
	#echo "j: "$j

#done

#base_dir="../results/tomtom" #TODO

#extract all motif names and paths
#Extract filename column from ../../TFBS/Results/tomtom/tomtom/metadata.tsv
