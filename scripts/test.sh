root=/fh/fast/setty_m/user/dotto/benchmarkDA
meld_bin=$root/methods/meld/bm_meld.py
cna_bin="$root/methods/cna/bm_cna.py"
data_id="bcr-xl"
data_id="covid19-pbmc"
data_dir=${root}/data/real/$data_id
pop=CD14_Monocyte
batch_sd=0
k=30
resolution=0.5
beta=25
downsample=3
mem=32g
seed=43
pop_enr=0.75
method=cydar

cmd="python $meld_bin \
--data_dir ${data_dir} \
--data_id ${data_id} \
--pop_enr $pop_enr \
--k $k \
--beta ${beta} \
--pop ${pop} \
--be_sd $batch_sd \
--seed $seed \
--outdir ${root}/benchmark/${data_id}/"

cmd="python $cna_bin \
--data_dir ${data_dir} \
--data_id ${data_id} \
--pop_enr $pop_enr \
--k $k \
--pop ${pop} \
--be_sd $batch_sd \
--seed $seed \
--outdir ${root}/benchmark/${data_id}/"
echo "$cmd"
