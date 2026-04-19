# See README.md for run instructions and expected outputs.
# See LICENSE.md and CONTRIBUTING.md for license and contribution details.
# this script clusters groups (k-means) of genome based on their distances (mash)
# and select representatives for each cluster
# also produce a plot of the clusters and the selected representatives
import logging
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.manifold import MDS

logger = logging.getLogger('scripts/cluster_and_select_genomes.py')
logging.basicConfig(level=logging.INFO)

# parameters are passed as positional arguments from Nextflow.
max_representatives = int(sys.argv[1])  # the number of clusters to create
group_name = sys.argv[2]
metadata_path = sys.argv[3]
use_ncbi = sys.argv[4] == 'true'
use_gtdb = sys.argv[5] == 'true'
fasta_files = sys.argv[6:]  # all the fasta to use as a list


def _write_annotated_fasta(
    input_fasta: str,
    output_fasta: str,
    header: str,
) -> None:
    """
    Write a new fasta file with the same sequence as the input fasta
    but with an annotated header.
    @param: input_fasta: the path to the input fasta file
    @param: output_fasta: the path to the output fasta file to write
    @param: header: the new header to write in the output fasta file
    """
    with open(output_fasta, 'w') as outp:
        for line in open(input_fasta):
            if line.startswith('>'):
                line = header
            outp.write(line)


def _get_taxids(metadata_path: str, group_name: str) -> tuple[str, str]:
    """
    Get the NCBI and GTDB taxids of a group of genomes from the metadata file.
    @param: metadata_path: the path to the metadata file
    @param: group_name: the name of the group of genomes to get the taxids for
    """
    taxid = '0'
    acc = 'NA'
    with open(metadata_path) as gtdb_metadata:
        for line in gtdb_metadata:
            if line.strip().split('\t')[-2] == group_name:
                taxid = line.split('\t')[79]
                acc = line.split('\t')[56]
                break

    gtdb_taxid = '0'
    with open('names.dmp') as gtdb_taxo_names:
        for line in gtdb_taxo_names:
            if acc in line:
                gtdb_taxid = line.split('\t')[0]
                break

    return taxid, gtdb_taxid


def _count_lines(path: str) -> int:
    """
    @param: path: the path to the file to count the lines of
    """
    n_lines = 0
    for _ in open(path):
        n_lines += 1
    return n_lines


def main() -> None:
    taxid, gtdb_taxid = _get_taxids(metadata_path, group_name)
    distance_tsv = f'{group_name}.distances.tsv'
    n_lines = _count_lines(distance_tsv)

    ncbi_state = pd.read_csv(
        metadata_path,
        usecols=['ncbi_genbank_assembly_accession', 'ncbi_keep'],
        sep='\t',
    )
    ncbi_state = ncbi_state[ncbi_state['ncbi_keep'] == 1]
    ncbi_keep_accessions = {
        value.split('.')[0]
        for value in ncbi_state[
            'ncbi_genbank_assembly_accession'
        ].values.tolist()
    }

    if n_lines < 3:
        for fasta in fasta_files:
            if use_ncbi and fasta.split('.')[0] in ncbi_keep_accessions:
                _write_annotated_fasta(
                    fasta,
                    f'{fasta[:-4]}.selected.ncbi.fna',
                    f'>{fasta[:-4]}|kraken:taxid|{taxid}\n',
                )
            if use_gtdb:
                _write_annotated_fasta(
                    fasta,
                    f'{fasta[:-4]}.selected.gtdb.fna',
                    f'>{fasta[:-4]}|kraken:taxid|{gtdb_taxid}|gtdb_taxonomy\n',
                )
        open(f'{group_name}.plot.png', 'w').close()
        raise SystemExit(0)

    df = pd.read_csv(
        distance_tsv,
        sep='\t',
        names=['D1', 'D2', 'distance', 'P', 'Cov'],
        dtype={'A': str, 'B': str, 'distance': np.float64},
    )
    logger.info('df loaded')

    df['D1'] = df['D1'].apply(lambda row: re.findall('GCA_[0-9]*', row)[0])
    df['D2'] = df['D2'].apply(lambda row: re.findall('GCA_[0-9]*', row)[0])
    logger.info('label shortened')

    logger.info('pivoting...')
    d_matrix = df.pivot(index='D1', columns='D2', values='distance')
    logger.info('pivoted')

    d_matrix = pd.concat(
        [
            pd.DataFrame(
                pd.Series([0.0], index=[d_matrix.iloc[:, 0].name]).add(
                    d_matrix.iloc[:, 0], fill_value=0
                ),
                columns=[d_matrix.iloc[:, 0].name],
            ).transpose(),
            d_matrix,
        ],
        sort=False,
    )
    d_matrix.iloc[0, :] = 0.0
    d_matrix.fillna(0.0, inplace=True)
    d_matrix = d_matrix + d_matrix.transpose()
    logger.info('equilibrated')

    logger.info('MDS...')
    embedding = MDS(n_components=2, dissimilarity='precomputed')
    transformed = embedding.fit_transform(d_matrix)
    transformed_df = pd.DataFrame(
        transformed,
        index=d_matrix.index,
        columns=['Dimension 1', 'Dimension 2'],
        dtype=float,
    )
    logger.info('done...')

    logger.info('Starting kmeans...')
    n_clusters = min(max_representatives, len(transformed))
    if len(transformed) < max_representatives:
        logger.info('reducing n to %s...', n_clusters)

    kmeans = KMeans(n_clusters=n_clusters, random_state=0)
    kmeans.fit(transformed)

    selected_points = []
    for centroid in kmeans.cluster_centers_:
        distances = np.linalg.norm(transformed - centroid, axis=1)
        closest_index = np.argmin(distances)
        selected_points.append(transformed[closest_index])

    for point in selected_points:
        for fasta in fasta_files:
            is_selected = str(
                transformed_df[
                    (transformed_df['Dimension 1'] == point[0])
                    & (transformed_df['Dimension 2'] == point[1])
                ].index[0]
            ) in fasta
            if not is_selected:
                continue

            if use_ncbi and fasta.split('.')[0] in ncbi_keep_accessions:
                _write_annotated_fasta(
                    fasta,
                    f'{fasta[:-4]}.selected.ncbi.fna',
                    f'>{fasta[:-4]}|kraken:taxid|{taxid}\n',
                )
            if use_gtdb:
                _write_annotated_fasta(
                    fasta,
                    f'{fasta[:-4]}.selected.gtdb.fna',
                    f'>{fasta[:-4]}|kraken:taxid|{gtdb_taxid}|gtdb_taxonomy\n',
                )

    plt.scatter(
        transformed[:, 0],
        transformed[:, 1],
        color='blue',
        label='All points',
    )
    plt.scatter(
        kmeans.cluster_centers_[:, 0],
        kmeans.cluster_centers_[:, 1],
        color='red',
        label='Selected points',
    )
    plt.legend()
    plt.savefig(f'{group_name}.plot.png')
    logger.info('All done...')


if __name__ == '__main__':
    main()
