#!/bin/bash
array=$1 #this varies between 0 and 1000


cd /projappl/project_2006203/MOSTA-SSTAT/Experiments
# cd csc_scratch/motifsimilarity-private/experiments/

results_path=/scratch/project_2006203/MOSTA-SSTAT
mkdir $results_path
mkdir $results_path"/results_final_version2.2_correct"

#n=3854
#n*(n-1)/2 #7424731

#n=4003
#n*(n-1)/2 #8010003

#n=3993
#n*(n-1)/2 #7970028

#n=3933 #CORRECT
#n*(n-1)/2 #7732278

#7420000/1000=7420  
#1000*7420 
#(1000+1)*7420

#8010000/1000=8010  
#1000*8010
#(1000+1)*8010

#7970000/1000 7970

#7730000/1000 #7730

start_ind=$(($array*7730))
end_ind=$((($array+1)*7730 ))

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
  filenames=( $(cut -d ',' -f1 ../../TFBS/PWMs_final_version2.2/filenames_similarity_computations.csv ) )
  length=${#filenames[@]}
  TF_PATH=../../TFBS/
  TF1=${filenames[$i]}
  TF2=${filenames[$j]}

  #convert "pwms" to "transfac"
  TF1="${TF1/pwms/transfac}"
  TF2="${TF2/pwms/transfac}"

  echo $TF1
  echo $TF2


  motifnames=( $(cut -d ',' -f1 ../../TFBS/PWMs_final_version2.2/motifnames.csv ) )

  PWM1=${motifnames[$i]}
  PWM2=${motifnames[$j]}

  echo "../../TFBS/"$TF1 > $LOCAL_SCRATCH"/"$PWM1"_"$PWM2".list" #overwrite
  echo "../../TFBS/"$TF2 >> $LOCAL_SCRATCH"/"$PWM1"_"$PWM2".list" #add

  if [[ $index -eq $start_ind ]]
  then
     #echo $index is equal to start
     
     ../sstat .5 list:$LOCAL_SCRATCH"/"$PWM1"_"$PWM2".list" typeI 0.01 > $LOCAL_SCRATCH"/result_"$array".out"



  else
     #echo $index is greater then start
     ../sstat .5 list:$LOCAL_SCRATCH"/"$PWM1"_"$PWM2".list" typeI 0.01 >> $LOCAL_SCRATCH"/result_"$array".out"


  fi

  rm $LOCAL_SCRATCH/$PWM1"_"$PWM2".list"


done

echo "Done"
cd $LOCAL_SCRATCH
ls
cp "result_"$array".out" $results_path"/results_final_version2.2_correct/"


