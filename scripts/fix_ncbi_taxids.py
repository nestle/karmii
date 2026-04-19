# See README.md for run instructions and expected outputs.
# See LICENSE.md and CONTRIBUTING.md for license and contribution details.
# this script takes the metadata file and the NCBI taxonomy dump files
# to fix the NCBI taxids in the metadata file

import sys

import pandas as pd

# parameters are passed as positional arguments from Nextflow.
metadata_path = sys.argv[1]
output_path = sys.argv[2]
dmp_path = sys.argv[3]  # it should be merged.dmp


def main() -> None:
    metadata = pd.read_csv(metadata_path, sep='\t', low_memory=False)

    with open(dmp_path) as handle:
        replacements = {
            int(line.strip().split('\t')[0]): int(line.strip().split('\t')[1])
            for line in handle
        }

    metadata['ncbi_species_taxid'] = metadata['ncbi_species_taxid'].replace(
        replacements
    )
    metadata['ncbi_taxid'] = metadata['ncbi_taxid'].replace(replacements)
    metadata.to_csv(output_path, sep='\t', index=False)


if __name__ == '__main__':
    main()
