#!/usr/bin/env python3

# See README.md for instructions on how to run the pipeline and the expected outputs.
# See LICENSE.md and CONTRIBUTING.md for details on the license and how to contribute to the project.

import sys
import glob
import matplotlib.pyplot as plt
import numpy as np
from concurrent.futures import ThreadPoolExecutor
import seaborn as sns
colorblind_palette = sns.color_palette("colorblind", int(sys.argv[1]))


def calculate_confidence(read):
    total_taxid = 0
    total_0 = 0
    parts = read.split()
    for part in parts:
        if "|" in part:
            continue
        if part.startswith("0:"):
            total_0 += int(part.split(':')[1])
        elif part.startswith("A:"):
            pass
        else:
            total_taxid += int(part.split(':')[1])
    if total_taxid <= 1:
        total_taxid = 0
    if total_0 == 0:
        ratio = 1
    try:
        ratio = total_taxid / (total_taxid+total_0)
    except ZeroDivisionError:
        ratio = 0
    return ratio


def parallel_function(species_reads, thresholds):
    print(f"taxonomic_confidence:loading {species_reads}...")
    try:
        reads = [line.strip().split(
            '\t')[4] for line in open(species_reads)]
        print(f"taxonomic_confidence:{species_reads} loaded, processing...")
    except IndexError:
        print(species_reads)
    read_confidences = {}
    for read in reads:
        read_confidences[read] = calculate_confidence(read)
    species_confidences = []
    with open(f"{species_reads}.conf.txt", 'w') as outp:
        for threshold in thresholds:
            states = []
            for read in reads:
                states.append('U' if (
                    read_confidences[read] < threshold/100
                    or read_confidences[read] == 0) else 'C')
            species_confidence = len(
                [x for x in states if x == 'C'])/len(
                    states)
            outp.write(f"{threshold/100:.4f},{species_confidence:.4f}\n")
            species_confidences.append(species_confidence)
            if species_reads not in max_confs:
                max_confs[species_reads] = species_confidence
        outp.write(
            f"# slope:"
            f"{((species_confidences[1]-species_confidences[10])/0.5)*100}\n")
        slopes[species_reads] = ((
            species_confidences[1]-species_confidences[10])/0.5)*100
    print(f"taxonomic_confidence:{species_reads} done")


thresholds = list(range(0, 100, 5))+[100]

print('taxonomic_confidence:Evaluating the confidence '
      'of Kraken2 species predictions...')
with ThreadPoolExecutor(max_workers=int(sys.argv[2])) as executor:
    max_confs = {}
    slopes = {}
    list(executor.map(
        lambda species_reads: parallel_function(species_reads, thresholds),
        glob.glob("*reads_classification*.tsv")))
print("taxonomic_confidence:parallel steps done...")

# produce a plot
# Create a figure
plt.figure()
fig = plt.figure()
ax = plt.subplot(111)

# Loop through each file and plot its data
n = 0
for i, species_reads in enumerate(
    sorted(max_confs, key=max_confs.get, reverse=True)):
    file = f"{species_reads}.conf.txt"
    for line in open(file):
        if line.startswith("#"):
            slope_value = float(line.split(":")[1])
    n += 1
    # Load the data from the text file
    data = np.loadtxt(file, delimiter=",", comments='#')
    x = data[:, 0]  # First column for x values
    y = data[:, 1]  # Second column for y values

    # Plot the data
    name = "N/A"
    for line in open(glob.glob(
        f"*.{species_reads.split('.')[-3]}.single_profile.tsv")[0]):
        if line.strip().split('\t')[3] == 'S':
            name = line.strip().split('\t')[5].lstrip()
            break
    if slopes[species_reads] > 20:
        alpha = 0.33
        linestyle = '--'
    else:
        alpha = 1
        linestyle = '-'
    if name != "N/A":
        plt.plot(x, y,
        label=f"{name} (-{slopes[species_reads]:.0f}%)",
        alpha=alpha, linestyle=linestyle,
        color=colorblind_palette[i])
if n == 0:
    print('taxonomic_confidence:ERROR:no file to process')

# Add titles and labels
plt.xlabel('Confidence')
plt.ylabel('Percent classified')
plt.ylim(-0.05, 1.05)
plt.legend()  # Show legend with file names

ax.axhline(y=0.95, color='lightgrey', linewidth=0.5, linestyle='-')
ax.axvline(x=0.05, color='lightgrey', linewidth=0.5, linestyle='-')
ax.axvline(x=0.5, color='lightgrey', linewidth=0.5, linestyle='-')

# Shrink current axis by 60%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.4, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.4), handletextpad=0.2)

# Save the plot as a PNG file
plt.savefig(f"{sys.argv[3]}.{sys.argv[4]}.confidence_plot.png", dpi=600)
