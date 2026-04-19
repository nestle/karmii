## Summary

This Nextflow pipeline analyzes microbial single‑isolate sequencing reads (FASTQ format; paired‑end or single‑end) using a two‑step workflow.

Its purpose is to assign a taxonomic label to the sample using either the NCBI or GTDB taxonomy and to detect potential contamination by secondary organisms at the read level. The workflow also generates a visual representation of the classification confidence for each candidate species.

The first step performs a standard Kraken analysis with a multi‑species database to generate an initial profile of the possible composition. The second step consists of a species‑specific analysis of the best candidate species using single‑species databases.

The pipeline includes a module for constructing and curating a custom reference database, which should be run once before the first analysis.

## License - please refer to LICENSE.md
This repository is licensed under the Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0).  
In summary, you may use, modify, and share this code for non-commercial purposes only, with proper attribution.  
Commercial use is prohibited except by copyright holder.

Full license text: [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/)

## Contributions - please refer to CONTRIBUTING.md
By contributing to this repository, you agree to the Contributor License Agreement (CLA), which grants the copyright holder of the original work in this repository full rights to use and commercialize your contributions.  

## External Dependencies
This pipeline uses third-party tools (e.g., kraken2, mash, matplotlib, pandas, scikit-learn, seaborn) under their respective licenses.  
Please consult their documentation for license terms. This repository does **not** redistribute these tools or any external data.

## Userguide
This code contains two Nextflow pipelines and thus requires Nextflow (https://nextflow.io/) to be installed and able to work with conda (https://anaconda.org/).

`nextflow run prepare_database.nf -resume` will collect all genomes necessary for building the reference database. Prior to this, you will need to:

1. Place one GTDB metadata file in `db_preparation_inputs/` following the pattern `gtdb_metadata_*.tsv`. This should contain all the organisms you want to have in your database (you can concat multiple files, for instance, if you want bacteria and archaea). For tests, try using a file limited to ~1000 lines.

> For instance https://data.gtdb.aau.ecogenomic.org/releases/release226/226.0/bac120_metadata_r226.tsv.gz (need to be unzipped and renamed gtdb_metadata_r226.tsv)
> 
2. Place a NCBI taxonomy dmp archive in `db_preparation_inputs/`following the pattern `*taxdump.tar.gz`.
> For instance https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz

3. Edit the `nextflow.config` file, which is documented to be self-explanatory.

Pay attention that this parameter remains set to `false` when you first run the `prepare_database.nf` workflow.

```yaml
params.downloads_validated = false
```

4. Proceed to a first run with `nextflow run prepare_database.nf -resume`

> Refer to `https://www.nextflow.io/docs/latest/` for details about how to customize your Nextflow run

4. (b) Some downloads may fail due to corrupted datasets or other reasons on the NCBI server side. Try rerunning `nextflow run prepare_database.nf -resume` a few times as some failures may not be systematic.

> Note: if you want to minimize the missing genomes, you can set `params.max_dl_per_request` to 1, although this will increase the time to download genomes dramatically.

4. (c) At some point, some downloads will still fail, set the following parameter to true to indicate you are done trying to obtain all genomes

```yaml
params.downloads_validated = true
```

5. Run `nextflow run prepare_database.nf -resume` again and let the database preparation process finish.

6. In the `krakendbs/ncbi/` and `krakendbs/gtdb/` folders, manually build the Kraken2 databases as follows:

```bash
kraken2-build --threads n --build --db .
```

> This is not automated with the nf worflow as it is rather resources consuming and careful evaluation of the resources is necessary for this step. When done, it will not affect the nf pipeline status so it will not rerun if nothing else has changed.

> If you make changes that will retrigger the `prepare_database.nf` pipeline, first remove the `krakendbs/ncbi/` and `krakendbs/gtdb/` folders to ensure the next kraken2 database will be built using only the files selected by the last run of the nf pipeline (published folders in nf are not cleaned during a run). In doubt, always remove these folders and rerun the pipeline, it will recreate from the cache a clean versions of these folders, ready for a new manual call to `kraken2-build`.

> If you do not have a Kraken2 install, you can load the conda environment created by the pipeline in the Nextflow work directory. Refer to the nf documentation to locate the work directory.

```bash
find /path/to/nf/work/conda -type f -name "kraken2"
```

**Once this is done, you are ready to run analyses**

7. Place your samples, either paired end (matching the pattern `samples/*_R{1,2}*.fastq.gz`) or single end (`samples/*.fastq.gz` excluding R1 and R2) in the samples subfolder.

8. Run 

```bash
nextflow run analysis.nf -resume
```

9. You will get your analysis results, as follows

```bash
samples/
├── mysample_R1_001.fastq.gz
└── mysample_R2_001.fastq.gz
results_profiling/
├── multi # results from first step, comprehensive db
│   ├── profile
│   │   ├── mysample.gtdb.multi_profile.tsv
│   │   └── mysample.ncbi.multi_profile.tsv
│   └── reads
│       ├── mysample.gtdb.reads_classification.tsv
│       └── mysample.ncbi.reads_classification.tsv
└── single # results from 2nd step, for each species, here n_species=2
    ├── profile
    │   ├── mysample.gtdb.s__Salmonella_enterica.single_profile.tsv
    │   ├── mysample.gtdb.s__Yersinia_pestis.single_profile.tsv
    │   ├── mysample.ncbi.28901.single_profile.tsv
    │   └── mysample.ncbi.632.single_profile.tsv
    └── reads
        ├── mysample.gtdb.s__Salmonella_enterica.reads_classification.tsv
        ├── mysample.gtdb.s__Yersinia_pestis.reads_classification.tsv
        ├── mysample.ncbi.28901.reads_classification.tsv
        └── mysample.ncbi.632.reads_classification.tsv
confidence_plots/ # plot summary of the analysis
├── mysample.gtdb.confidence_plot.png
└── mysample.ncbi.confidence_plot.png
├── mysample.gtdb.confidence_plot.svg
└── mysample.ncbi.confidence_plot.svg
```

**Example of confidence plot for one sample:**

![Example confidence plot for a single sample. The plot shows 10 lines, with the y‑axis representing the percentage of classified reads and the x‑axis representing the Kraken confidence threshold. Listeria monocytogenes shows a prediction starting at 95% and maintains this value as the confidence threshold increases. The other tested species start at 60% or lower and drop immediately as the confidence threshold increases.](images/confidence_plot.png)

> Signature of a clean isolate of *B. longum*. The percentage reported for each line represents the slope of the curve between confidence x=0.05 and confidence x=0.5.