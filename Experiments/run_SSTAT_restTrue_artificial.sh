#!/bin/bash
array=$1 #this varies between 0 and 1000


cd /projappl/project_2006203/MOSTA-SSTAT/Experiments
# cd csc_scratch/motifsimilarity-private/experiments/

results_path=/scratch/project_2006203/MOSTA-SSTAT/results_restTrue_artificial
mkdir $results_path

#1031*24744=25511064
#1031*24744/1000=25511.06

num_of_rows=1031
num_of_columns=24744

#Test different number of cpus:
#start_ind=$(($array*10)) 
#end_ind=$((($array+1)*10 ))

start_ind=$(($array*25511)) 
end_ind=$((($array+1)*25511 ))

if [ "$end_ind" -gt 25511063 ]; then
   
   end_ind=25511063
fi

readarray -t representatives < /projappl/project_2006203/TFBS/PWMs_final/representatives.csv

filenames=( $(cut -d ',' -f1 ../../TFBS/PWMs_final/filenames.csv ) ) #These are the full filenames of all motifs
TF_PATH=../../TFBS/

#True motif
readarray -t filenames_artificial < /projappl/project_2006203/TFBS/PWMs_final/representatives_artificial.csv

artificial_motif_path="/scratch/project_2006203/TFBS/artificial_motifs/artificial_motifs_transfac/"


for (( index=$start_ind; index<$end_ind; index++ )); do

    row=$((index / num_of_columns))
    col=$((index % num_of_columns))
    
    true_motif=${representatives[$row]}
    true_motif_file=$(printf "%s\n" "${filenames[@]}" | grep "$true_motif")
    
    #echo ${#filenames[@]} 24744
    
    
    TF1=$true_motif_file
    #convert "pwms" to "transfac"
    TF1="${TF1/pwms/transfac}"
    
    TF2=${filenames_artificial[$col]}
    
    F1="${TF1/pwms/transfac}"
    
    PWM1=$(basename "$TF1")
    PWM1="${PWM1%.*}"
    
    PWM2=$(basename "$TF2")
    TF2=$artificial_motif_path$PWM2
    PWM2="${PWM2%.*}"
    
    echo "../../TFBS/"$TF1 > $LOCAL_SCRATCH"/"$PWM1"_"$PWM2".list"
    echo $TF2 >> $LOCAL_SCRATCH"/"$PWM1"_"$PWM2".list"
    
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

cd $LOCAL_SCRATCH
cp "result_"$array".out" $results_path"/"


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
