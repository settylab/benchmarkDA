#!/bin/bash

# add slurm module first
module purge
module load R/4.3.1-gfbf-2022b
module load ImageMagick/7.1.0-53-GCCcore-12.2.0
module load GSL/2.7-GCCcore-12.2.0
module load cuDNN/8.4.1.50-CUDA-11.7.0
eval "$(micromamba shell hook --shell=bash)"

# set slurm parameters
time=1-00:00:00
partition=campus-new

data_id=$1
script_path="$(readlink -f "$0")"
script_dir="$(dirname "$script_path")"
if [ -z "${root+x}" ]; then
    export root="$(readlink -f "$script_dir/..")"
fi
cd ${root}/scripts

if [[ "$data_id" == "cluster" ]]
    then
    data_dir=${root}/data/synthetic/$data_id
    pops=$(for p in $(seq 1 1 3); do echo M$p; done)
    R_methods=$(for m in mellon mellon_dm milo daseq cydar cna meld louvain milo_batch cna_batch louvain_batch; do echo $m; done)
    batch_vec=$(for m in 0 0.75 1 1.25 1.5; do echo $m; done)
    k=30
    resolution=0.2
    beta=33
    downsample=3
    mem=8g
elif [[ "$data_id" == "cluster_balanced" ]]
    then
    data_dir=${root}/data/synthetic/$data_id
    pops=$(for p in $(seq 1 1 3); do echo M$p; done)
    R_methods=$(for m in mellon mellon_dm milo daseq cydar cna meld louvain milo_batch cna_batch louvain_batch; do echo $m; done)
    batch_vec=$(for m in 0 0.75 1 1.25 1.5; do echo $m; done)
    k=30
    resolution=0.2
    beta=33
    downsample=3
    mem=8g
elif [[ "$data_id" == "linear" ]]
    then
    data_dir=${root}/data/synthetic/$data_id
    pops=$(for p in $(seq 1 1 7); do echo M$p; done)
    R_methods=$(for m in mellon mellon_dm milo daseq cydar cna meld louvain milo_batch cna_batch louvain_batch; do echo $m; done)
    batch_vec=$(for m in 0 0.75 1 1.25 1.5; do echo $m; done)
    k=30
    resolution=1
    beta=71
    downsample=3
    mem=8g
elif [[ "$data_id" == "branch" ]]
    then
    data_dir=${root}/data/synthetic/$data_id
    pops=$(for p in $(seq 1 1 8); do echo M$p; done)
    R_methods=$(for m in mellon mellon_dm milo daseq cydar cna meld louvain milo_batch cna_batch louvain_batch; do echo $m; done)
    batch_vec=$(for m in 0 0.75 1 1.25 1.5; do echo $m; done)
    k=30
    resolution=1
    beta=65
    downsample=3
    mem=8g
elif [[ "$data_id" == "covid19-pbmc" ]]
    then
    data_dir=${root}/data/real/$data_id
    pops=$(for m in RBC B PB CD14_Monocyte CD8_T CD4_T Platelet NK Granulocyte CD16_Monocyte gd_T pDC DC; do echo $m; done)
    R_methods=$(for m in mellon mellon_dm milo daseq cydar cna meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolution=0.5
    beta=25
    downsample=3
    mem=32g
elif [[ "$data_id" == "bcr-xl" ]]
    then
    data_dir=${root}/data/real/$data_id
    pops=$(for m in CD4_T-cells NK_cells CD8_T-cells B-cells_IgM+ monocytes surface- B-cells_IgM- DC; do echo $m; done)
    R_methods=$(for m in mellon mellon_dm milo daseq cydar cna meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolution=0.6
    beta=23
    downsample=10
    mem=32g
fi


job_number=0
## Run
for pop in $pops
    do
    for pop_enr in $(seq 0.75 0.1 0.95)
        do
        for seed in 43 44 45
            do
            for batch_sd in $batch_vec
                do
                for method in $R_methods
                    do
                    ((job_number++))
                    if [ -z "$SLURM_ARRAY_TASK_ID" ] || [ "$job_number" -ne "$SLURM_ARRAY_TASK_ID" ]; then
                        continue
                    fi
                    jobid=${data_id}-${pop}-${pop_enr}-${batch_sd}-${method}-${seed}
                    echo "Doing $jobid ..."
                    if [[ "$method" == "cna" ]]; then
                        # enalbe cna env
                        micromamba activate cna
                        cna_bin="$root/methods/cna/bm_cna.py"
                        python $cna_bin \
                            --data_dir ${data_dir} \
                            --data_id ${data_id} \
                            --pop_enr $pop_enr \
                            --k $k \
                            --pop ${pop} \
                            --be_sd $batch_sd \
                            --seed $seed \
                            --outdir ${root}/benchmark/${data_id}/
                        exit $!
                    elif [[ "$method" == "cna_batch" ]]; then
                        # enalbe cna env
                        micromamba activate cna
                        cna_bin="$root/methods/cna/bm_cna.py"
                        python $cna_bin \
                            --data_dir ${data_dir} \
                            --data_id ${data_id} \
                            --pop_enr $pop_enr \
                            --k $k \
                            --pop ${pop} \
                            --be_sd $batch_sd \
                            --seed $seed \
                            --outdir ${root}/benchmark/${data_id}/ \
                            --model_batch
                        exit $!
                    elif [[ "$method" == "meld" ]]; then
                        # enable meld env
                        micromamba activate meld
                        meld_bin="$root/methods/meld/bm_meld.py"
                          python $meld_bin \
                            --data_dir ${data_dir} \
                            --data_id ${data_id} \
                            --pop_enr $pop_enr \
                            --k $k \
                            --beta ${beta} \
                            --pop ${pop} \
                            --be_sd $batch_sd \
                            --seed $seed \
                            --outdir ${root}/benchmark/${data_id}/
                        exit $!
                    elif [[ "$method" == "mellon" ]]; then
                        # enable meld env
                        micromamba activate mellon_v2
                        meld_bin="$root/methods/mellon/bm_mellon.py"
                          python $meld_bin \
                            --data_dir ${data_dir} \
                            --data_id ${data_id} \
                            --pop_enr $pop_enr \
                            --pop ${pop} \
                            --dm-comp 0 \
                            --be_sd $batch_sd \
                            --seed $seed \
                            --out-name $method \
                            --outdir ${root}/benchmark/${data_id}/
                        exit $!
                    elif [[ "$method" == "mellon_dm" ]]; then
                        # enable meld env
                        micromamba activate mellon_v2
                        meld_bin="$root/methods/mellon/bm_mellon.py"
                          python $meld_bin \
                            --ls-factor 100 \
                            --data_dir ${data_dir} \
                            --data_id ${data_id} \
                            --pop_enr $pop_enr \
                            --pop ${pop} \
                            --dm-comp 10 \
                            --be_sd $batch_sd \
                            --seed $seed \
                            --out-name $method \
                            --outdir ${root}/benchmark/${data_id}/
                        exit $!
                    elif [[ "$method" == "mellon_hls" ]]; then
                        # enable meld env
                        micromamba activate mellon_v2
                        meld_bin="$root/methods/mellon/bm_mellon.py"
                          python $meld_bin \
                            --ls-factor 100 \
                            --data_dir ${data_dir} \
                            --data_id ${data_id} \
                            --pop_enr $pop_enr \
                            --pop ${pop} \
                            --be_sd $batch_sd \
                            --seed $seed \
                            --out-name $method \
                            --outdir ${root}/benchmark/${data_id}/
                        exit $!
                    else
                        Rscript ./run_DA.r \
                            ${data_dir}/${data_id}_data_bm.RDS $method $seed $pop \
                            --data_dir ${data_dir}/ \
                            --pop_enrichment $pop_enr \
                            --data_id $data_id \
                            --k $k \
                            --resolution ${resolution} \
                            --downsample ${downsample} \
                            --batchEffect_sd $batch_sd \
                            --outdir ${root}/benchmark/${data_id}/
                        exit $!
                    fi
                done
            done
        done
    done
done

# Submit a slurm array job
jobid="benchmarkDA_syn_real_$data_id"
cmd="sbatch -J '$jobid' --time=$time --partition=$partition \
--mem $mem --out '$root/SlurmLog/${jobid}_%N_%A_%a.out' --array=1-$job_number \
'$script_path' $1"
echo "$cmd"
eval "$cmd"
