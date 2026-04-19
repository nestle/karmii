/*
 * See README.md for instructions on how to run the pipeline and the expected outputs.
 * See LICENSE.md and CONTRIBUTING.md for license and contribution details.
 */

 /*
 * This nf pipeline prepares the databases for Kraken2 profiling.
 * It takes as input the GTDB metadata and the NCBI taxonomy dump,
 * filters the metadata according to user-defined quality thresholds,
 * downloads the corresponding genomes from NCBI,
 * clusters them by group,
 * selects representatives from each cluster,
 * and adds them to the Kraken2 databases.
 */

/*
 * Channels
 */
Channel
    .fromPath('db_preparation_inputs/gtdb_metadata_*.tsv')
    .set { raw_metadata }

Channel
    .fromPath('db_preparation_inputs/*taxdump.tar.gz')
    .set { taxdump }

Channel
    .fromPath('scripts/metadata_filtering.py')
    .set { metadata_filtering_script }

Channel
    .fromPath('scripts/ncbi_bad_states_filtering.py')
    .set { ncbi_bad_states_filtering_script }

Channel
    .fromPath('scripts/cluster_and_select_genomes.py')
    .set { cluster_and_select_genomes_script }

Channel
    .fromPath('scripts/fix_ncbi_taxids.py')
    .set { fix_ncbi_taxids_script }

params.downloaded = "${workflow.workDir}/downloaded.txt"

/*
 * Processes
 */
process metadata_filtering {
    label 'big_task'
    publishDir 'metadata/', mode: params.publish_mode, overwrite: true
    conda 'conda/metadata_filtering.yaml'

    input:
    path raw_metadata
    path pythonScript

    output:
    path 'metadata_filtered.tsv'

    script:
    """
    python3 ${pythonScript} ${raw_metadata} metadata_filtered.tsv \
        ${params.checkm_minimum_completeness} \
        ${params.checkm_maximum_contamination} \
        ${params.exclude_envs_and_mag} \
        ${params.max_nb_genomes_per_group_to_process} 2>&1
    """
}

process get_genome_ncbi_states {
    maxForks params.apikey ? params.apikeynjobs : 1
    publishDir 'metadata/ncbi_states/', mode: params.publish_mode, overwrite: true
    conda 'conda/get_genome_ncbi_states.yaml'

    input:
    val filtered_accessions

    output:
    path "${filtered_accessions[0]}.json"

    script:
    """
    datasets summary genome accession \
        ${filtered_accessions.join(' ').replace('GCA_', 'GCF_')} \
        ${filtered_accessions.join(' ')} \
        ${params.apikey ? '--api-key ' + params.apikey : ''} \
        > ${filtered_accessions[0]}.json
    """
}

process ncbi_bad_states_filtering {
    label 'big_task'
    publishDir 'metadata/', mode: params.publish_mode, overwrite: true
    conda 'conda/ncbi_bad_states_filtering.yaml'

    input:
    path filtered_accessions
    path genome_ncbi_states
    path pythonScript

    output:
    path 'metadata_ncbi_state_annotated.tsv'

    script:
    """
    python3 ${pythonScript} ${filtered_accessions} metadata_ncbi_state_annotated.tsv \
        ${genome_ncbi_states} 2>&1
    """
}

process download_genomes {
    errorStrategy 'ignore'
    maxForks params.apikey ? params.apikeynjobs : 1
    conda 'conda/download_genomes.yaml'

    input:
    val groups_and_accessions

    output:
    tuple path('*.fna'), val("${groups_and_accessions[0]}")

    when:
    !file(params.downloaded).exists() ||
        params.downloads_validated == false ||
        file(params.downloaded).text.trim().split('\n')
            .contains(groups_and_accessions[1..-1].join('_'))

    script:
    """
    datasets download genome accession ${groups_and_accessions[1..-1].join(' ')} \
        ${params.apikey ? '--api-key ' + params.apikey : ''} --no-progressbar --dehydrated
    unzip -n ncbi_dataset.zip
    datasets rehydrate --directory .
    mv ncbi_dataset/data/*/*.fna .
    chmod 660 *.fna
    rm ncbi_dataset.zip
    rm -r ncbi_dataset
    echo ${groups_and_accessions[1..-1].join('_')} >> ${params.downloaded}
    """
}

process mash_sketch {
    conda 'conda/mash_sketch.yaml'

    input:
    path(genomes)
    val(groups)

    output:
    tuple path('*.fna.msh'), val(groups)

    script:
    """
    mash sketch ${genomes}
    """
}

process mash_triangle {
    label 'big_task'
    conda 'conda/mash_triangle.yaml'

    input:
    tuple path(sketches), val(groups)

    output:
    tuple path('*.distances.tsv'), val(groups)

    script:
    """
    mash triangle -p ${task.cpus} -E *.msh > ${groups}.distances.tsv
    """
}

process prepare_ncbi_taxonomy {
    publishDir 'krakendbs/ncbi/taxonomy/', mode: params.publish_mode, overwrite: true

    input:
    path taxdump

    output:
    path '*.dmp'

    script:
    """
    tar zxvf ${taxdump}
    """
}

process prepare_gtdb_taxonomy {
    publishDir 'krakendbs/gtdb/taxonomy/', mode: params.publish_mode, overwrite: true
    conda 'conda/prepare_gtdb_taxonomy.yaml'

    input:
    path raw_metadata

    output:
    tuple path('*.dmp'), path('*.tsv')

    script:
    """
    python -c "import pandas as pd; df = pd.read_csv('${raw_metadata}', sep='\t'); df[['ncbi_genbank_assembly_accession', 'gtdb_taxonomy']].to_csv('metadata.tmp', sep='\t', index=False, header=False)"
    gtdb_to_taxdump.py metadata.tmp > taxID_info.tsv
    """
}

process fix_ncbi_taxonomy {
    publishDir 'metadata/', mode: params.publish_mode, overwrite: true
    conda 'conda/fix_ncbi_taxonomy.yaml'

    input:
    path dmp
    path metadata_filtered_twice
    path pythonScript

    output:
    path 'metadata_final.tsv'

    script:
    """
    cut -f1,3 merged.dmp > taxid_to_replace.tsv
    python3 ${pythonScript} ${metadata_filtered_twice} metadata_final.tsv taxid_to_replace.tsv
    """
}

process cluster_and_select_genomes {
    label 'big_task'
    publishDir "clusters/${groups}", mode: params.publish_mode, overwrite: true
    conda 'conda/cluster_and_select_genomes.yaml'

    input:
    tuple (
        val(groups),
        path(genomes),
        path(distances),
        path(pythonScript),
        path(metadata_ready),
        path(gtdb_dmp)
    )

    output:
    tuple path('*.selected.ncbi.fna'), path('*.selected.gtdb.fna'), path('*.plot.png')

    script:
    """
    python3 ${pythonScript} ${params.max_nb_representatives_per_group_to_keep} ${groups} \
        ${metadata_ready} ${params.ncbi} ${params.gtdb} ${genomes.collect().join(' ')} 2>&1
    touch none.selected.ncbi.fna
    touch none.selected.gtdb.fna
    """
}

process krakendb_ncbi_add {
    errorStrategy 'ignore'
    publishDir 'krakendbs/ncbi/', mode: params.publish_mode, overwrite: true
    conda 'conda/kraken.yaml'

    input:
    path(genomes)

    output:
    path 'library/added/*.txt'
    path 'library/added/*.masked'
    path 'library/added/*.fna'

    script:
    """
    kraken2-build --add-to-library ${genomes} --db .
    """
}

process krakendb_gtdb_add {
    errorStrategy 'ignore'
    publishDir 'krakendbs/gtdb/', mode: params.publish_mode, overwrite: true
    conda 'conda/kraken.yaml'

    input:
    path(genomes)

    output:
    path 'library/added/*.txt'
    path 'library/added/*.masked'
    path 'library/added/*.fna'

    script:
    """
    kraken2-build --add-to-library ${genomes} --db .
    """
}

/*
 * Workflow
 */
workflow {
    if (params.ncbi || params.gtdb) {
        metadata_filtering(raw_metadata.first(), metadata_filtering_script)
            .set { metadata_filtered_once }

        metadata_filtered_once
            .splitCsv(header: true, sep: '\t')
            .map { row -> row.'ncbi_genbank_assembly_accession' }
            .set { filtered_accessions }

        get_genome_ncbi_states(
            filtered_accessions.buffer(
                size: (params.max_json_per_query).toInteger(),
                remainder: true
            )
        ).set { genomes_ncbi_states }

        // mark entries that are not taxonomy OK in NCBI
        ncbi_bad_states_filtering(
            metadata_filtered_once,
            genomes_ncbi_states.collect(),
            ncbi_bad_states_filtering_script
        ).set { metadata_filtered_twice }

        metadata_filtered_twice
            .splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.'group', row.'ncbi_genbank_assembly_accession') }
            .set { groups_and_ncbi_genbank_assembly_accessions }

        // group accessions by group and split into batches of max_dl_per_request
        // to avoid overwhelming NCBI servers.
        groups_and_ncbi_genbank_assembly_accessions
            .groupTuple()
            .map { batch -> tuple(batch[0], batch[1].collate(params.max_dl_per_request)) }
            .set { ncbi_genbank_assembly_accessions_by_group }

        // repeat group numbers according to the number of batches of accessions
        // to download for each group
        ncbi_genbank_assembly_accessions_by_group
            .map { item -> [item[0]] * item[1].size() }
            .flatten()
            .set { repeated_group_numbers }

        // flatten batches of accessions to download for each group
        ncbi_genbank_assembly_accessions_by_group
            .map { item -> item[1] }
            .flatMap { it }
            .set { accessions }

        // merge repeated group numbers with accessions to get a tuple of (group, accession)
        repeated_group_numbers
            .merge(accessions)
            .set { groups_and_accessions }

        download_genomes(groups_and_accessions)
            .set { downloaded_genomes }

        if (params.downloads_validated) {
            // when only one genome is downloaded for a group, the output is a tuple of (path, group)
            // when multiple genomes are downloaded for a group, the output is a tuple of (list of paths, group)
            // this conditional ensures that the output is always a tuple of (list of paths, group) to simplify downstream processing.
            downloaded_genomes
                .map { item ->
                    if (item[0].getClass().getName() == 'sun.nio.fs.UnixPath') {
                        return [[item[0]], item[1]]
                    } else {
                        return [item[0], item[1]]
                    }
                }
                .set { genomes_and_groups }

            // get rid of the multiple levels caused by max_dl_per_request
            // and reconstitute the original grouping of genomes by group
            genomes_and_groups
                .groupTuple(by: 1)
                .map { [it[0].flatten(), it[1]] }
                .set { genomes_regrouped }

            genomes_regrouped
                .map { it[0] }
                .set { genomes }

            genomes_regrouped
                .map { [it[1]] * it[0].size() }
                .set { groups }

            mash_sketch(genomes.flatten(), groups.flatten())
                .set { sketches_and_groups }

            sketches_and_groups
                .groupTuple(by: 1)
                .set { sketches_regrouped }

            mash_triangle(sketches_regrouped)
                .set { distances_grouped }

            prepare_ncbi_taxonomy(taxdump)
                .set { ncbi_dmp }

            // fix NCBI taxonomy to take into account merged.dmp
            // and get the final metadata file
            fix_ncbi_taxonomy(ncbi_dmp, metadata_filtered_twice, fix_ncbi_taxids_script)
                .set { metadata_ready }

            prepare_gtdb_taxonomy(metadata_ready)
                .map { item -> item[0] }
                .set { gtdb_dmp }

            cluster_and_select_genomes(
                genomes_regrouped
                    .combine(distances_grouped, by: 1)
                    .combine(cluster_and_select_genomes_script)
                    .combine(metadata_ready)
                    .combine(gtdb_dmp.flatten().filter { file -> file =~ /.*(names\.dmp)/ })
            ).set { clustered_and_selected_genomes }

            clustered_and_selected_genomes
                .map { it[0] }
                .flatten()
                .set { genomes_for_kraken_ncbi }

            clustered_and_selected_genomes
                .map { it[1] }
                .flatten()
                .set { genomes_for_kraken_gtdb }

            if (params.ncbi) {
                krakendb_ncbi_add(
                    genomes_for_kraken_ncbi.filter { file -> !file.name.equals('none.selected.ncbi.fna') }
                )
            }

            if (params.gtdb) {
                krakendb_gtdb_add(
                    genomes_for_kraken_gtdb.filter { file -> !file.name.equals('none.selected.gtdb.fna') }
                )
            }
        }
    }
}
