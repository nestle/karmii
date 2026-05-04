**v1.0.3**
- Filter genomic files matching unwanted patterns (e.g. _cds_from_genomic, _rna_from_genomic)
- Use `--no-masking` when building species specific dbs to ensure no k-mer is discarded
- [debug] Reorganise the workflow to avoid rerunning all single dbs analyses when a new db is produced

**v1.0.2**
- Use the ignore errorStrategy for the get_genome_ncbi_states process

**v1.0.1**
- Added an SVG version of the plot.
- Enabled the use of MAGs in the GTDB mode (can be disabled in the configuration).
- Use dehydrated mode for downloads, which should reduce the number of failures.

**v1.0.0**
- Initial release.