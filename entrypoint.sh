#!/bin/bash
source /opt/conda/etc/profile.d/conda.sh
conda activate snakemake
cd /usr/src/vcf/
exec "$@"