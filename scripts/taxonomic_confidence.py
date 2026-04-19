#!/usr/bin/env python3

# See README.md for run instructions and expected outputs.
# See LICENSE.md and CONTRIBUTING.md for license and contribution details.
# this script explore the kraken reads
# classification and produce the confidence plot

import glob
import sys
import textwrap
from concurrent.futures import ThreadPoolExecutor

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# parameters are passed as positional arguments from Nextflow.
nb_species = sys.argv[1]
cpus = sys.argv[2]
sample = sys.argv[3]
taxonomy = sys.argv[4]

thresholds = list(range(0, 100, 5)) + [100]
colorblind_palette = sns.color_palette('colorblind', int(nb_species))


def calculate_confidence(read: str) -> float:
    """
    Calculate the confidence above which a read is classifed
    @param: read: the string containing the kraken2 classification of a read
    """
    total_taxid = 0
    total_0 = 0

    for part in read.split():
        if '|' in part:
            continue
        if part.startswith('0:'):
            total_0 += int(part.split(':')[1])
        elif part.startswith('A:'):
            pass
        else:
            total_taxid += int(part.split(':')[1])

    if total_taxid <= 1:
        total_taxid = 0

    if total_0 == 0:
        ratio = 1
    else:
        ratio = 0

    try:
        ratio = total_taxid / (total_taxid + total_0)
    except ZeroDivisionError:
        ratio = 0

    return ratio


def species_confidence(
    species_reads: str,
    thresholds: list[int],
    max_confs: dict[str, float],
    slopes: dict[str, float],
) -> None:
    """
    Process a single species reads classification file to calculate confidence
    @param: species_reads: the path to the kraken2
    classification file of a species
    @param: thresholds: the list of confidence thresholds to evaluate
    @param: max_confs: a dictionary to store the maximum
    confidence of each species (used to sort species in the final plot)
    @param: slopes: a dictionary to store the slope of
    the confidence curve of each species
    """
    print(f'taxonomic_confidence:loading {species_reads}...')
    try:
        reads = [line.strip().split('\t')[4] for line in open(species_reads)]
        print(f'taxonomic_confidence:{species_reads} loaded, processing...')
    except IndexError:
        print(species_reads)
        reads = []

    read_confidences = {read: calculate_confidence(read) for read in reads}
    species_confidences = []

    with open(f'{species_reads}.conf.txt', 'w') as outp:
        for threshold in thresholds:
            states = [
                'U'
                if (
                    read_confidences[read] < threshold / 100
                    or read_confidences[read] == 0
                )
                else 'C'
                for read in reads
            ]
            species_confidence = (
                len([state for state in states if state == 'C']) / len(states)
            )
            outp.write(f'{threshold / 100:.4f},{species_confidence:.4f}\n')
            species_confidences.append(species_confidence)
            if species_reads not in max_confs:
                max_confs[species_reads] = species_confidence

        slope = (
            (species_confidences[1] - species_confidences[10]) / 0.5
        ) * 100
        outp.write(f'# slope:{slope}\n')
        slopes[species_reads] = slope

    print(f'taxonomic_confidence:{species_reads} done')


def wrap_label(label: str, width: int = 40) -> str:
    """
    Wrap a label to a specified width for better readability in plots.
    @param: label: the label string to wrap
    @param: width: the maximum width of each line
    """
    return textwrap.fill(label, width=width, break_on_hyphens=False)


def main() -> None:
    print(
        'taxonomic_confidence:Evaluating the confidence '
        'of Kraken2 species predictions...'
    )

    max_confs: dict[str, float] = {}
    slopes: dict[str, float] = {}
    with ThreadPoolExecutor(max_workers=int(cpus)) as executor:
        list(
            executor.map(
                lambda species_reads: species_confidence(
                    species_reads,
                    thresholds,
                    max_confs,
                    slopes,
                ),
                glob.glob('*reads_classification*.tsv'),
            )
        )

    print('taxonomic_confidence:parallel steps done...')

    plt.figure()
    ax = plt.subplot(111)

    n = 0
    for i, species_reads in enumerate(
        sorted(max_confs, key=max_confs.get, reverse=True)
    ):
        file_path = f'{species_reads}.conf.txt'
        for line in open(file_path):
            if line.startswith('#'):
                slope_value = float(line.split(':')[1])
                _ = slope_value

        n += 1
        data = np.loadtxt(file_path, delimiter=',', comments='#')
        x = data[:, 0]
        y = data[:, 1]

        name = 'N/A'
        profile_files = glob.glob(
            f'*.{species_reads.split(".")[-3]}.single_profile.tsv'
        )
        for line in open(profile_files[0]):
            if line.strip().split('\t')[3] == 'S':
                name = line.strip().split('\t')[5].lstrip()
                break

        if slopes[species_reads] > 20:
            alpha = 0.33
            linestyle = '--'
        else:
            alpha = 1
            linestyle = '-'

        if name != 'N/A':
            wrapped_label = wrap_label(
                f'{name} (-{slopes[species_reads]:.0f}%)',
                width=40,
            )
            plt.plot(
                x,
                y,
                label=wrapped_label,
                alpha=alpha,
                linestyle=linestyle,
                color=colorblind_palette[i],
            )

    if n == 0:
        print('taxonomic_confidence:ERROR:no file to process')

    plt.xlabel('Confidence')
    plt.ylabel('Percent classified')
    plt.ylim(-0.05, 1.05)
    plt.legend()

    ax.axhline(y=0.95, color='lightgrey', linewidth=0.5, linestyle='-')
    ax.axvline(x=0.05, color='lightgrey', linewidth=0.5, linestyle='-')
    ax.axvline(x=0.5, color='lightgrey', linewidth=0.5, linestyle='-')

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.4, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.4), handletextpad=0.2)

    plt.savefig(f'{sample}.{taxonomy}.confidence_plot.png', dpi=600)
    plt.savefig(f'{sample}.{taxonomy}.confidence_plot.svg')


if __name__ == '__main__':
    main()
