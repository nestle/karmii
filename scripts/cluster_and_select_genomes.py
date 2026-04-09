# See README.md for instructions on how to run the pipeline and the expected outputs.
# See LICENSE.md and CONTRIBUTING.md for details on the license and how to contribute to the project.

import pandas as pd
import numpy as np
import re
from sklearn.manifold import MDS
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

import sys
import logging

logger = logging.getLogger("scripts/cluster_and_select_genomes.py")
logging.basicConfig(level=logging.INFO)

# get taxid corresponding to the group
taxid = '0'
acc = 'NA'
with open(sys.argv[3]) as gtdb_metadata:
    for line in gtdb_metadata:
        if line.strip().split('\t')[-2] == sys.argv[2]:
            taxid = line.split('\t')[79]
            acc = line.split('\t')[56]
            break

with open('names.dmp') as gtdb_taxo_names:
    for line in gtdb_taxo_names:
        if acc in line:
            gtdb_taxid = line.split('\t')[0]

nl = 0
for line in open("{}.distances.tsv".format(sys.argv[2])):
    nl += 1

# load the ncbi_state of the genomes
ncbi_state = pd.read_csv(sys.argv[3],
                         usecols=['ncbi_genbank_assembly_accession',
                                  'ncbi_keep'], sep="\t")
ncbi_state = ncbi_state[ncbi_state['ncbi_keep'] == 1]
if nl < 3:
    for f in sys.argv[6:]:
        if sys.argv[4] == "true":
            if f.split('.')[0] in [v.split('.')[0] for v in ncbi_state[
                    'ncbi_genbank_assembly_accession'].values.tolist()]:
                with open('{}.selected.ncbi.fna'.format(f[:-4]), 'w') as outp:
                    for line in open(f):
                        if line.startswith(">"):
                            line = '>{}|kraken:taxid|{}\n'.format(
                                f[:-4], taxid)
                        outp.write(line)
        if sys.argv[5] == "true":
            with open('{}.selected.gtdb.fna'.format(f[:-4]), 'w') as outp:
                for line in open(f):
                    if line.startswith(">"):
                        line = f">{f[:-4]}|" \
                               f"kraken:taxid|{gtdb_taxid}|gtdb_taxonomy\n"
                    outp.write(line)
    open('{}.plot.png'.format(sys.argv[2]), 'w').close()
    exit(0)

df = pd.read_csv("{}.distances.tsv".format(sys.argv[2]), sep='\t',
                 names=['D1', 'D2', 'distance', 'P', 'Cov'],
                 dtype={'A': str, 'B': str, 'distance': np.float64})
logger.info('df loaded')
df["D1"] = df["D1"].apply(lambda row: re.findall('GCA_[0-9]*', row)[0])
df["D2"] = df["D2"].apply(lambda row: re.findall('GCA_[0-9]*', row)[0])
logger.info('label shortened')
logger.info('pivoting...')
# Pivot the DataFrame into a matrix
d_matrix = df.pivot(index='D1', columns='D2', values='distance')
logger.info('pivoted')

# fix what need to be fixed
# transpose the first column as first line
d_matrix = pd.concat([pd.DataFrame(pd.Series(
    [0.0], index=[d_matrix.iloc[:, 0].name]).add(
    d_matrix.iloc[:, 0], fill_value=0),
    columns=[d_matrix.iloc[:, 0].name]).transpose(), d_matrix], sort=False)
# but remove values
d_matrix.iloc[0, :] = 0.0
d_matrix.fillna(0.0, inplace=True)
# make an equilibrated matrix
logger.info('equilibrating')
d_matrix = d_matrix+d_matrix.transpose()
logger.info('equilibrated')

# Perform MDS
logger.info('MDS...')
embedding = MDS(n_components=2, dissimilarity='precomputed')
X_transformed = embedding.fit_transform(d_matrix)
df_transformed = pd.DataFrame(
    X_transformed, index=d_matrix.index,
    columns=['Dimension 1', 'Dimension 2'], dtype=float)
logger.info('done...')

logger.info('Starting kmeans...')
n = int(sys.argv[1])
if len(X_transformed) < n:
    n = len(X_transformed)
    logger.info('reducing n to {}...'.format(n))
kmeans = KMeans(n_clusters=n, random_state=0)

# Fit the KMeans algorithm to the transformed points
kmeans.fit(X_transformed)

# Get the cluster centroids
selected_indices = kmeans.labels_

selected_points = []

# Loop through each centroid
for centroid in kmeans.cluster_centers_:
    # Calculate the distances from each data point to the centroid
    distances = np.linalg.norm(X_transformed - centroid, axis=1)

    # Find the index of the data point
    # with the minimum distance to the centroid
    closest_index = np.argmin(distances)

    # Get the closest data point
    closest_point = X_transformed[closest_index]

    # Append the closest data point to the list
    selected_points.append(closest_point)

for p in selected_points:
    for f in sys.argv[6:]:
        if str(df_transformed[(df_transformed['Dimension 1'] == p[0])
                              & (df_transformed['Dimension 2'] == p[1])
                              ].index[0]) in f:
            if sys.argv[4] == "true":
                if f.split('.')[0] in [v.split('.')[0] for v in ncbi_state[
                         'ncbi_genbank_assembly_accession'].values.tolist()]:
                    with open('{}.selected.ncbi.fna'.format(
                            f[:-4]), 'w') as outp:
                        for line in open(f):
                            if line.startswith(">"):
                                line = '>{}|kraken:taxid|{}\n'.format(
                                    f[:-4], taxid)
                            outp.write(line)
            if sys.argv[5] == "true":
                with open('{}.selected.gtdb.fna'.format(f[:-4]), 'w') as outp:
                    for line in open(f):
                        if line.startswith(">"):
                            line = f">{f[:-4]}" \
                                   f"|kraken:taxid|" \
                                   f"{gtdb_taxid}|gtdb_taxonomy\n"
                        outp.write(line)

# plot the graph

# Assume X_transformed are your points in 2D after MDS
# and selected_points are the points you selected

# Create a scatter plot of all points
plt.scatter(X_transformed[:, 0], X_transformed[:, 1], color='blue',
            label='All points')

# Create a scatter plot of the selected points
plt.scatter(kmeans.cluster_centers_[:, 0], kmeans.cluster_centers_[:, 1],
            color='red', label='Selected points')

# Add a legend
plt.legend()

# Save the plot as a PNG file
plt.savefig('{}.plot.png'.format(sys.argv[2]))
logging.info('All done...')
