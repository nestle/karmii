#!/usr/bin/env python3

# See README.md for run instructions and expected outputs.
# See LICENSE.md and CONTRIBUTING.md for license and contribution details.
# This script filters the GTDB metadata file to mark genomes
# that do not pass the NCBI taxonomy checks

import json
import sys

import pandas as pd

# parameters are passed as positional arguments from Nextflow.
metadata_path = sys.argv[1]
output_path = sys.argv[2]
ncbi_states_json_paths = sys.argv[3:]


def main() -> None:
    metadata = pd.read_csv(metadata_path, sep='\t', low_memory=False)
    to_drop = set()

    taxids_dict = {
        row['ncbi_genbank_assembly_accession']: [
            int(row['ncbi_species_taxid']),
            int(row['ncbi_taxid']),
        ]
        for _, row in metadata.iterrows()
    }

    for jfile_name in ncbi_states_json_paths:
        with open(jfile_name, 'r') as jfile:
            jdata = json.load(jfile)
            for genome in jdata['reports']:
                #  drop if not one these two states
                if genome['assembly_info']['assembly_status'] not in [
                    'current',
                    'previous',
                ]:
                    to_drop.add(genome['accession'].replace('GCF_', 'GCA_'))

                try:
                    #  drop if taxo check towards type strain failed
                    #  or if this check is not available
                    if (
                        genome['average_nucleotide_identity'][
                            'taxonomy_check_status'
                        ]
                        not in ['OK']
                    ):
                        to_drop.add(
                            genome['accession'].replace('GCF_', 'GCA_')
                        )
                except KeyError:
                    to_drop.add(genome['accession'].replace('GCF_', 'GCA_'))

                try:
                    # drop if the taxid is not present
                    if int(genome['organism']['tax_id']) not in taxids_dict[
                        genome['accession']
                    ]:
                        to_drop.add(
                            genome['accession'].replace('GCF_', 'GCA_')
                        )
                except KeyError:
                    accession = genome['accession'].replace('GCF_', 'GCA_')
                    if int(genome['organism']['tax_id']) not in taxids_dict[
                        accession
                    ]:
                        to_drop.add(accession)

    metadata['ncbi_keep'] = 1
    metadata.loc[
        metadata['ncbi_genbank_assembly_accession'].isin(to_drop),
        'ncbi_keep',
    ] = 0
    metadata.to_csv(output_path, sep='\t', index=False)


if __name__ == '__main__':
    main()
