#!/bin/bash
array=$1 #this varies between 0 and 393


cd /projappl/project_2006203/MOSTA-SSTAT/Experiments
# cd csc_scratch/motifsimilarity-private/experiments/

results_path=/scratch/project_2006203/MOSTA-SSTAT
mkdir $results_path
mkdir $results_path"/results_true_vs_corresponding_artificial_version2.2"

readarray -t motifs < /projappl/project_2006203/TFBS/PWMs_final_version2.2/motifnames.csv

#filenames=( $(cut -d ',' -f1 ../../TFBS/PWMs_final/filenames.csv ) ) #These are the full filenames of all motifs

filenames=( $(cut -d ',' -f1 ../../TFBS/PWMs_final_version2.2/filenames_similarity_computations.csv ) ) #These are the full filenames of all motifs


start_ind=$(($array*10)) 
end_ind=$((($array+1)*10 ))

if [ "$end_ind" -gt 3933 ]; then
   
   end_ind=3933
fi

write_new_file=1
search_dir="/scratch/project_2006203/TFBS/PWMs_final_version2.2/artificial_motifs_transfac"

for (( index=$start_ind; index<$end_ind; index++ )); do

 true_motif=${motifs[$index]}".pfm"
 
 true_motif_file=$(printf "%s\n" "${filenames[@]}" | grep "$true_motif")
  
 

  # Empty array to hold the resulting files
  artificial_files=()

  true_motif=${motifs[$index]}

  # For each prefix/lines use glob to get files and add them to the result array
  for file in "$search_dir/$true_motif"*; do
      # Check if file exists (glob might return pattern if no files are found)
      if [[ -f $file ]]; then
          artificial_files+=("$file")
      fi
  done
  
  
   if [[ "$true_motif" == *"v2"* ]]; then
   #echo "The string '$true_motif' contains the substring v2. Do nothing"
    :
   else
    #echo "The string '$true_motif' does not contain the substring v2. Remove possible v2-motifs from the artificial list"
    substring="v2"
    # Print each element of the array on a new line, filter with grep -v (invert match), and read into a new array
    mapfile -t filtered_artificial_files < <(printf "%s\n" "${artificial_files[@]}" | grep -v "$substring")
    # Print the filtered_array to verify
    #printf "%s\n" "${filtered_artificial_files[@]}"
    artificial_files=("${filtered_artificial_files[@]}")
  fi
  
  
  
  PWM1=$true_motif
  
  TF1=$true_motif_file

  #convert "pwms" to "transfac"
  TF1="${TF1/pwms/transfac}"
  
  
  for TF2 in "${artificial_files[@]}"; do #Loop over artificial
  
    PWM2=$(basename "$TF2")
    PWM2="${PWM2%.*}"
    
   
    echo "../../TFBS/"$TF1 > $LOCAL_SCRATCH"/"$PWM1"_"$PWM2".list"
    echo $TF2 >> $LOCAL_SCRATCH"/"$PWM1"_"$PWM2".list"
    
    if [[ $write_new_file -eq 1 ]]
    then
       #echo $index is equal to start
      
    
       ../sstat .5 list:$LOCAL_SCRATCH"/"$PWM1"_"$PWM2".list" typeI 0.01 > $LOCAL_SCRATCH"/result_"$array".out"
       write_new_file=0


    else
       #echo $index is greater then start
       
       ../sstat .5 list:$LOCAL_SCRATCH"/"$PWM1"_"$PWM2".list" typeI 0.01 >> $LOCAL_SCRATCH"/result_"$array".out"
    fi

    rm $LOCAL_SCRATCH/$PWM1"_"$PWM2".list"
    
  
  done #Loop over artificial
  
done #Loop over index
  
cd $LOCAL_SCRATCH
cp "result_"$array".out" $results_path"/results_true_vs_corresponding_artificial_version2.2/"

