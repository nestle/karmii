# See README.md for run instructions and expected outputs.
# See LICENSE.md and CONTRIBUTING.md for license and contribution details.
# this script filters the GTDB metadata file
# to keep only genomes matching our quality criteria

import logging
import random
import string
import sys

import pandas as pd

logger = logging.getLogger('scripts/metadata_filtering.py')
logging.basicConfig(level=logging.INFO)

# parameters are passed as positional arguments from Nextflow.
metadata_path = sys.argv[1]
output_path = sys.argv[2]
completeness_threshold = sys.argv[3]
contamination_threshold = sys.argv[4]
ignore_mag_and_env = sys.argv[5]
max_to_process_per_group = sys.argv[6]

# fix columns to gtdb r226 to ensure the rest of the nf script based
# on cut -f position of columns works as expected
# if columns are added in the future.
columns = ['accession',
           'ambiguous_bases',
           'checkm2_completeness',
           'checkm2_contamination',
           'checkm2_model',
           'checkm_completeness',
           'checkm_contamination',
           'checkm_marker_count',
           'checkm_marker_lineage',
           'checkm_marker_set_count',
           'checkm_strain_heterogeneity',
           'coding_bases',
           'coding_density',
           'contig_count',
           'gc_count',
           'gc_percentage',
           'genome_size',
           'gtdb_genome_representative',
           'gtdb_representative',
           'gtdb_taxonomy',
           'gtdb_type_designation_ncbi_taxa',
           'gtdb_type_designation_ncbi_taxa_sources',
           'gtdb_type_species_of_genus',
           'l50_contigs',
           'l50_scaffolds',
           'longest_contig',
           'longest_scaffold',
           'lsu_23s_contig_len',
           'lsu_23s_count',
           'lsu_23s_length',
           'lsu_23s_query_id',
           'lsu_5s_contig_len',
           'lsu_5s_count',
           'lsu_5s_length',
           'lsu_5s_query_id',
           'lsu_silva_23s_blast_align_len',
           'lsu_silva_23s_blast_bitscore',
           'lsu_silva_23s_blast_evalue',
           'lsu_silva_23s_blast_perc_identity',
           'lsu_silva_23s_blast_subject_id',
           'lsu_silva_23s_taxonomy',
           'mean_contig_length',
           'mean_scaffold_length',
           'mimag_high_quality',
           'mimag_low_quality',
           'mimag_medium_quality',
           'n50_contigs',
           'n50_scaffolds',
           'ncbi_assembly_level',
           'ncbi_assembly_name',
           'ncbi_assembly_type',
           'ncbi_bioproject',
           'ncbi_biosample',
           'ncbi_contig_count',
           'ncbi_contig_n50',
           'ncbi_country',
           'ncbi_date',
           'ncbi_genbank_assembly_accession',
           'ncbi_genome_category',
           'ncbi_genome_representation',
           'ncbi_isolate',
           'ncbi_isolation_source',
           'ncbi_lat_lon',
           'ncbi_molecule_count',
           'ncbi_ncrna_count',
           'ncbi_organism_name',
           'ncbi_protein_count',
           'ncbi_refseq_category',
           'ncbi_rrna_count',
           'ncbi_scaffold_count',
           'ncbi_scaffold_l50',
           'ncbi_scaffold_n50',
           'ncbi_scaffold_n75',
           'ncbi_scaffold_n90',
           'ncbi_seq_rel_date',
           'ncbi_spanned_gaps',
           'ncbi_species_taxid',
           'ncbi_ssu_count',
           'ncbi_strain_identifiers',
           'ncbi_submitter',
           'ncbi_taxid',
           'ncbi_taxonomy',
           'ncbi_taxonomy_unfiltered',
           'ncbi_total_gap_length',
           'ncbi_total_length',
           'ncbi_translation_table',
           'ncbi_trna_count',
           'ncbi_type_material_designation',
           'ncbi_ungapped_length',
           'ncbi_unspanned_gaps',
           'ncbi_wgs_master',
           'protein_count',
           'scaffold_count',
           'ssu_contig_len',
           'ssu_count',
           'ssu_gg_blast_align_len',
           'ssu_gg_blast_bitscore',
           'ssu_gg_blast_evalue',
           'ssu_gg_blast_perc_identity',
           'ssu_gg_blast_subject_id',
           'ssu_gg_taxonomy',
           'ssu_length',
           'ssu_query_id',
           'ssu_silva_blast_align_len',
           'ssu_silva_blast_bitscore',
           'ssu_silva_blast_evalue',
           'ssu_silva_blast_perc_identity',
           'ssu_silva_blast_subject_id',
           'ssu_silva_taxonomy',
           'total_gap_length',
           'trna_aa_count',
           'trna_count',
           'trna_selenocysteine_count']


def main() -> None:
    gtdb_md = pd.read_csv(metadata_path, sep='\t', low_memory=False)[columns]

    # drop columns that contains taxonomic annotations
    # from other sources than GTDB and NCBI
    # and might create confusion
    columns_to_drop = [
        'lsu_silva_23s_taxonomy',
        'ncbi_taxonomy_unfiltered',
        'ssu_gg_taxonomy',
        'ssu_silva_taxonomy',
    ]
    gtdb_md = gtdb_md.drop(columns=columns_to_drop)
    logger.info(f'Nb of genomes is {len(gtdb_md)}')

    if ignore_mag_and_env:
        gtdb_md = gtdb_md[
            ~gtdb_md.astype(str)
            .apply(lambda x: x.str.contains('derived from metagenome'))
            .any(axis=1)
        ]
        logger.info(
            'Nb of genomes is %s after removing derived from metagenome',
            len(gtdb_md),
        )

        gtdb_md = gtdb_md[
            ~gtdb_md.astype(str)
            .apply(lambda x: x.str.contains('-MAG'))
            .any(axis=1)
        ]
        logger.info(f'Nb of genomes is {len(gtdb_md)} after removing -MAG')

        gtdb_md = gtdb_md[
            ~gtdb_md.astype(str)
            .apply(lambda x: x.str.contains('uncultured'))
            .any(axis=1)
        ]
        logger.info(f'Nb of genomes is {len(gtdb_md)} '
                    f'after removing uncultured')

        gtdb_md = gtdb_md[~gtdb_md['ncbi_taxonomy'].str.endswith('s__')]
        logger.info(f'Nb of genomes is {len(gtdb_md)} after removing s__')

    gtdb_md = gtdb_md[
        gtdb_md['checkm_completeness'].astype(float) >= float(
            completeness_threshold)
    ]
    gtdb_md = gtdb_md[
        gtdb_md['checkm_contamination'].astype(float) <= float(
            contamination_threshold)
    ]
    logger.info(
        (
            'Nb of genomes is %s after removing checkM completeness < %s%% '
            'and contamination > %s%%'
        ),
        len(gtdb_md),
        completeness_threshold,
        contamination_threshold,
    )

    gtdb_md['group'] = (
        gtdb_md.groupby(['ncbi_taxonomy', 'gtdb_taxonomy']).ngroup() + 1
    )
    gtdb_md['group'] = gtdb_md['group'].astype(str) + '-' + ''.join(
        random.choices(string.ascii_lowercase, k=10)
    )
    logger.info(
        'Nb of group is %s once grouped by ncbi_taxonomy+gtdb_taxonomy',
        gtdb_md['group'].max(),
    )

    gtdb_md = gtdb_md.sort_values('group')
    gtdb_md = (
        gtdb_md.groupby('group')
        .apply(
            lambda x: x.sample(
                len(x) if len(x) < int(max_to_process_per_group)
                else int(max_to_process_per_group)
            )
        )
        .reset_index(drop=True)
    )

    gtdb_md.to_csv(output_path, sep='\t', index=False)
    logger.info('%s DONE %s', '-' * 10, '-' * 10)


if __name__ == '__main__':
    main()
