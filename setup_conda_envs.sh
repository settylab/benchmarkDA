# !/bin/bash -

micromamba -y -n meld -c conda-forge -c bioconda scanpy pip
micromamba activate meld
pip install meld
micromamba deactivate
create -n cna scanpy pip python-can
