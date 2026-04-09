# See README.md for instructions on how to run the pipeline and the expected outputs.
# See LICENSE.md and CONTRIBUTING.md for details on the license and how to contribute to the project.

import pandas as pd
import sys
import logging
import json

metadata = pd.read_csv(sys.argv[1], sep="\t", low_memory=False)

to_drop = set()

#
# Create a dictionnary to map the accession and taxids
taxids_dict = {row['ncbi_genbank_assembly_accession']:
               [int(row['ncbi_species_taxid']), int(row['ncbi_taxid'])]
               for _, row in metadata.iterrows()}
for jfile_name in sys.argv[3:]:
    with open(jfile_name, 'r') as jfile:
        jdata = json.load(jfile)
        for genome in jdata["reports"]:
            # check status
            if genome["assembly_info"]["assembly_status"] not in ['current',
                                                                  'previous']:
                to_drop.add(genome["accession"].replace('GCF_', 'GCA_'))

            # check the taxonomic validation done with ANI to type strain
            try:
                if genome["average_nucleotide_identity"][
                          "taxonomy_check_status"] not in ['OK']:
                    to_drop.add(genome["accession"].replace('GCF_', 'GCA_'))
            except KeyError:
                to_drop.add(genome["accession"].replace('GCF_', 'GCA_'))
            # the GTDB file is always behind the live NCBI state,
            # if the taxid has changed in NCBI compared to the initial file,
            # the metadata is wrong, so drop it
            try:
                if int(genome["organism"]
                       ["tax_id"]) not in taxids_dict[genome["accession"]]:
                    to_drop.add(genome["accession"].replace('GCF_', 'GCA_'))
            except KeyError:
                if int(genome["organism"]
                       ["tax_id"]) not in taxids_dict[
                           genome["accession"].replace('GCF_', 'GCA_')]:
                    to_drop.add(genome["accession"].replace('GCF_', 'GCA_'))

# mark accessions that are problematic.
# They are not dropped but marked since gtdb can still use them
metadata["ncbi_keep"] = 1
metadata.loc[metadata[
    'ncbi_genbank_assembly_accession'].isin(to_drop), "ncbi_keep"] = 0

# write
metadata.to_csv(sys.argv[2], sep="\t", index=False)
