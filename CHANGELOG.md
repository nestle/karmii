**v1.0.5**
- Add an installation section to the README
- Add a Dockerfile and adapt some parts of the code to be compatible with containerisation
- Add a summary to facilitate decision-making
**v1.0.4**
- Make the code compliant with strict Nextflow vocabulary (nf 26+)
- Support as input files *.fastq, *.fasta, *.fasta.gz in addition to *.fastq.gz
- Change of parameter (and logic) `n_species` to `n_genus_species`
- Fix a metadata filtering issue with recent pandas versions
- Ensure forced redownload when necessary without manual action
- Provide empty default input folders

**v1.0.3**
- Filter genomic files matching unwanted patterns (e.g. _cds_from_genomic, _rna_from_genomic)
- Use `--no-masking` when building species-specific DBs to ensure no k-mer is discarded
- [debug] Reorganise the workflow to avoid rerunning all single-DB analyses when a new db is produced

**v1.0.2**
- Use the ignore errorStrategy for the get_genome_ncbi_states process

**v1.0.1**
- Add an SVG version of the plot.
- Enable the use of MAGs in GTDB mode (can be disabled in the configuration).
- Use dehydrated mode for downloads, which should reduce the number of failures.

**v1.0.0**
- Initial release.