# See README.md for instructions on how to run the pipeline and the expected outputs.
# See LICENSE.md and CONTRIBUTING.md for details on the license and how to contribute to the project.

import pandas as pd
import sys
import logging
import random
import string

logger = logging.getLogger("scripts/gtdb_filtering.py")
logging.basicConfig(level=logging.INFO)
gtdb_md = pd.read_csv(sys.argv[1], sep="\t", low_memory=False)
# drop taxonomies that I don't want to use for filtering
columns_to_drop=["lsu_silva_23s_taxonomy", "ncbi_taxonomy_unfiltered", "ssu_gg_taxonomy", "ssu_silva_taxonomy"]
gtdb_md = gtdb_md.drop(columns=columns_to_drop)
logger.info(f"Nb of genomes is {len(gtdb_md)}")
gtdb_md = gtdb_md[~gtdb_md.astype(str).apply(lambda x: x.str.contains('derived from metagenome')).any(axis=1)]
logger.info(f"Nb of genomes is {len(gtdb_md)} after removing derived from metagenome")
gtdb_md = gtdb_md[~gtdb_md.astype(str).apply(lambda x: x.str.contains('-MAG')).any(axis=1)]
logger.info(f"Nb of genomes is {len(gtdb_md)} after removing -MAG")
gtdb_md = gtdb_md[~gtdb_md.astype(str).apply(lambda x: x.str.contains('uncultured')).any(axis=1)]
logger.info(f"Nb of genomes is {len(gtdb_md)} after removing uncultured")
gtdb_md = gtdb_md[~gtdb_md['ncbi_taxonomy'].str.endswith('s__')]
logger.info(f"Nb of genomes is {len(gtdb_md)} after removing s__")
gtdb_md = gtdb_md[gtdb_md['checkm_completeness'].astype(float) >= float(sys.argv[3])]
gtdb_md = gtdb_md[gtdb_md['checkm_contamination'].astype(float) <= float(sys.argv[4])]
logger.info(f"Nb of genomes is {len(gtdb_md)} after removing checkM completeness < {sys.argv[3]}% and contamination > {sys.argv[4]}%")
gtdb_md['group'] = gtdb_md.groupby(['ncbi_taxonomy', 'gtdb_taxonomy']).ngroup() + 1
gtdb_md['group'] = gtdb_md['group'].astype(str) + "-" + ''.join(random.choices(string.ascii_lowercase, k=10))
logger.info(f"Nb of group is {gtdb_md['group'].max()} once grouped by ncbi_taxonomy+gtdb_taxonomy")
# sort by group so all similar genome will be downloaded together if this order is used
gtdb_md = gtdb_md.sort_values('group')
# randomly select max x genomes for each groups
gtdb_md = gtdb_md.groupby('group').apply(lambda x: x.sample(len(x) if len(x) < int(sys.argv[5]) else int(sys.argv[5]))).reset_index(drop=True)
# write
gtdb_md.to_csv(sys.argv[2], sep="\t", index=False)
logger.info("-"*10+" DONE "+"-"*10)
