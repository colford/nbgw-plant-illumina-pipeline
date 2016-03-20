#!/bin/bash

# CD to the fasta/DNA directory
cd ../fasta/DNA
ls > ../../jobs/fasta-files.txt
cd -

# Read the fasta-files.txt
while read fasta
do
    ID=$fasta
    parts=(${ID//./ })
    cp blast-file-template.txt slurm/${parts[0]}.slurm
    sed_fasta="s/replaceme1/${ID}/"
    sed_csv="s/replaceme2/${parts[0]}/"
    perl -pi -e ${sed_fasta} slurm/${parts[0]}.slurm
    perl -pi -e ${sed_csv} slurm/${parts[0]}.slurm 
done < fasta-files.txt

# Keep count of nodes required
nodes=0

# Create the bsub submitter file
rm slurm/jobs_*.sh
numjobs=0
for fname in slurm/*.slurm 
do
    if [ $((nodes % 12)) -eq 0 ]
    then
    	numjobs=$((numjobs + 1)) 
	echo '#!/bin/bash' > slurm/jobs_${numjobs}.sh
    fi

    echo "sbatch ${fname}" >> slurm/jobs_${numjobs}.sh
    nodes=$((nodes + 1))
done
chmod +x slurm/jobs_*.sh

# Let me know how many nodes are required
echo "${nodes} nodes required to run all jobs"

# Now run the jobs in batches
for job_file in slurm/job_*.sh 
do
    echo "Running: ${job_file}"
    source ${job_file}
    while true 
    do
       str=`squeue|grep cpb128`
       if [[ -z $str ]]
       then
          break
       fi
       sleep 10
    done
done
