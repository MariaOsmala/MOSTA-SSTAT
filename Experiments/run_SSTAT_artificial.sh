#!/bin/bash
array=$1 #this varies between 0 and 1000


cd /projappl/project_2006203/MOSTA-SSTAT/Experiments
# cd csc_scratch/motifsimilarity-private/experiments/

results_path=/scratch/project_2006203/MOSTA-SSTAT
mkdir $results_path
mkdir $results_path"/artificial_results_final"
#FOXO1 FOXK1 esimerkki 624389

#start_ind=$(($array*8000))
#end_ind=$((($array+1)*8000 ))

#start_ind=$(($array*5420))
#end_ind=$((($array+1)*5420 ))

start_ind=$(($array*306121)) 
end_ind=$((($array+1)*306121 ))

readarray -t filenames < /projappl/project_2006203/TFBS/PWMs_final/representatives_artificial.csv

transfac_path=/scratch/project_2006203/TFBS/artificial_motifs/artificial_motifs_transfac/

for (( index=$start_ind; index<$end_ind; index++ )); do

  tmp1=$((1+8*$index))
  tmp=`echo sqrt "($tmp1)" | bc -l`
  tmp1=`echo "($tmp)+1" | bc -l`
  i=`echo "($tmp1)/2" | bc -l`
  i=${i%.*}
  j=$(($index-$i*($i-1)/2))

  #echo "index: "$index
  echo "i: "$i
  echo "j: "$j

#From i j back to index
#index=j+i*(i-1)/2


#filenames=( $(cut -d ',' -f1 ../../TFBS/Results/tomtom/tomtom/filenames.csv ) )
  
  #length=${#filenames[@]}
  #TF_PATH=../../TFBS/
  TF1=${filenames[$i]}
  TF2=${filenames[$j]}
  
  PWM1=$(basename "$TF1")
  TF1=$transfac_path$PWM1
  PWM1="${PWM1%.*}"
    
  PWM2=$(basename "$TF2")
  TF2=$transfac_path$PWM2
  PWM2="${PWM2%.*}"

  


#PWM1=${TF1#*${TF1%/*}}
#PWM1=${PWM1#*/}
#PWM1=${PWM1%.*}

#PWM2=${TF2#*${TF2%/*}}
#PWM2=${PWM2#*/}
#PWM2=${PWM2%.*}



  echo $TF1 > $LOCAL_SCRATCH"/"$PWM1"_"$PWM2".list"
  echo $TF2 >> $LOCAL_SCRATCH"/"$PWM1"_"$PWM2".list"

  if [[ $index -eq $start_ind ]]
  then
     echo $index is equal to start
     #echo $PWM1"-"$PWM2 > "../results/result_"$array".out"

     ../sstat .5 list:$LOCAL_SCRATCH"/"$PWM1"_"$PWM2".list" typeI 0.01 > $LOCAL_SCRATCH"/result_"$array".out"



  else
     echo $index is greater then start
     #echo $PWM1"-"$PWM2 >> "../results/result_"$array".out"
     ../sstat .5 list:$LOCAL_SCRATCH"/"$PWM1"_"$PWM2".list" typeI 0.01 >> $LOCAL_SCRATCH"/result_"$array".out"


  fi

  rm $LOCAL_SCRATCH/$PWM1"_"$PWM2".list"


done

cd $LOCAL_SCRATCH
cp "result_"$array".out" $results_path"/artificial_results_final/"


