from pathlib import Path
import warnings
import time
import argparse
import numpy as np
import pandas as pd
import os.path as osp
import logging

from _cna import runCNA, cna2output
from util import csv2anndata, calculate_outcome, out2bm, c_index

logging.basicConfig(
        level=logging.INFO, 
        format='[%(asctime)s] [%(levelname)-8s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
)


def parse_args():
    parser = argparse.ArgumentParser("benchmark the performance of cna method.")

    parser.add_argument("--data_dir", type=str, default="./data", help="data directory")
    parser.add_argument("--data_id", type=str, default="linear", help="the dataset id which specifies topology and reps")
    parser.add_argument("--k", type=int, default=30, help="K parameter to build KNN graph")
    parser.add_argument("--pop", type=str, default="M1", help="the population number")
    parser.add_argument("--pop_enr", type=str, default=0.75, help="population enrichment")
    parser.add_argument("--be_sd", type=float, default=0, help="SD of batch effect")
    parser.add_argument("--seed", type=int, default=43, help="random seed")
    parser.add_argument("--out_type", type=str, default="label", help="output data format")
    parser.add_argument("--outdir", type=str, help="dir to save benchmark results")
    parser.add_argument("--model_batch", action="store_true", help="whether add batch in the model")   

    args = parser.parse_args()
    return args

def runDA(args):
    # set random seed for reproducibility
    np.random.seed(args.seed)
    # read the data
    data_dir = args.data_dir
    prefix = "benchmark_{}_pop_{}_enr{}_seed{}".format(args.data_id,
                                                        args.pop,
                                                        args.pop_enr,
                                                        args.seed)
    pcaData = osp.join(data_dir, prefix + "_batchEffect{}.pca.csv".format(
        int(args.be_sd) if args.be_sd.is_integer() else args.be_sd))
    colData = osp.join(data_dir, prefix + ".coldata.csv")
    adata = csv2anndata(data_path=pcaData, obs_path=colData)
    
    # run analysis using CNA method
    logging.info("stating analysis")
    start_time = time.time()
    cna_res, md = runCNA(adata,
                         k=args.k,
                         sample_col="synth_samples",
                         label_col="synth_labels",
                         batch_col="synth_batches" if args.model_batch else None)
    run_time = time.time() - start_time
    
    cindex = c_index(cna_res.ncorrs, adata.obs['Condition2_prob'].values)
    
    if cna_res.p > .05:
        logging.warning("Global association p-value: {} > .05".format(cna_res.p))
    # get da cells with different alphas
    logging.info("preparing the result")
    alphas = np.percentile(cna_res.fdrs['fdr'], np.arange(1e-8, 1-1e-8, 0.01) * 100)
    da_cell = cna2output(md, cna_res, out_type=args.out_type, alphas=alphas) # out: [np.array([str])]
    
    # prepare benchmark data frame
    if args.out_type == "continuous":
        bm_out = out2bm(da_cell, adata, args.model_batch, run_time)
    else:
        bm_out = []
        for i, da in enumerate(da_cell):
            bm_out.append(calculate_outcome(out2bm(da, adata, args.model_batch)))
            bm_out[i][['method', 'alpha', 'runtime']] = ['cna', alphas[i], run_time]
            
    if args.out_type != "continuous":
        bm_out = pd.concat(bm_out, axis=0)
        # print("AUC: ", metrics.auc(bm_out['FPR'].values, bm_out['TPR'].values))
    
    # save the benchmark result
    outdir = args.outdir
    Path(outdir).mkdir(parents=True, exist_ok=True)
    cindex_file = osp.join(outdir, prefix + "_batchEffect{}.DAresults.{}".format(
        int(args.be_sd) if args.be_sd.is_integer() else args.be_sd, 'cna_batch' if args.model_batch else 'cna') + ".cindex")
    logging.info(f'writing c-index of {cindex} to "{cindex_file}".')
    with open(cindex_file, 'w') as f:
        f.write(str(cindex))
    logging.info("writing result csv")
    bm_resfile = osp.join(outdir, prefix + "_batchEffect{}.DAresults.{}".format(
        int(args.be_sd) if args.be_sd.is_integer() else args.be_sd, 'cna_batch' if args.model_batch else 'cna') + ".csv")
    bm_out.to_csv(bm_resfile, index=False)
    logging.info("success")


if __name__ == "__main__":
    args = parse_args()
    runDA(args)
