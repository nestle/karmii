# See README.md for instructions on how to run the pipeline and the expected outputs.
# See LICENSE.md and CONTRIBUTING.md for details on the license and how to contribute to the project.

import pandas as pd
import sys
md = pd.read_csv(sys.argv[1], sep="\t", low_memory=False)
replacements = {int(l.strip().split('\t')[0]):int(l.strip().split('\t')[1]) for l in open(sys.argv[3])}
md['ncbi_species_taxid'] = md['ncbi_species_taxid'].replace(replacements)
md['ncbi_taxid'] = md['ncbi_taxid'].replace(replacements)
# write
md.to_csv(sys.argv[2], sep="\t", index=False)