/*
* See README.md for instructions on how to run the pipeline and the expected outputs.
* See LICENSE.md and CONTRIBUTING.md for details on the license and how to contribute to the project.
*/

Channel
    .fromFilePairs("samples/*_R{1,2}*.fastq.gz", flat: true)
    .set { paired_samples }

Channel
    .fromPath("samples/*.fastq.gz")
    .filter { file -> 
        !(file.name =~ /.*_R1.*\.fastq.gz$/ || file.name =~ /.*_R2.*\.fastq.gz$/)
    }
    .set { single_samples }

Channel.fromPath("krakendbs/ncbi/")
    .set{ ncbi_multi_species_db }

Channel.fromPath("krakendbs/gtdb/")
    .set{ gtdb_multi_species_db }

Channel.fromPath("clusters/")
    .set{ genome_clusters }

Channel.fromPath("metadata/metadata_final.tsv")
    .set{ metadata }

Channel.fromPath("scripts/taxonomic_confidence.py")
    .set{ taxonomic_confidence_script }

process multi_profiling {
    label 'big_task'
    publishDir "results_profiling/multi/reads/", mode: params.publish_mode, overwrite: true, pattern: '*.reads_classification.tsv'
    publishDir "results_profiling/multi/profile/", mode: params.publish_mode, overwrite: true, pattern: '*.multi_profile.tsv'
    conda 'conda/kraken.yaml'
    input:
    tuple val(sample_id), path(r1), path(r2)
    each path(multi_species_db)
    output:
    tuple path("${sample_id}.${multi_species_db}.multi_profile.tsv"), val(multi_species_db), val(sample_id), path(r1), path(r2), path("${sample_id}.${multi_species_db}.reads_classification.tsv")
    script:
    
    if (r2.toString() == 'null')
    """
    kraken2 --threads ${task.cpus} --confidence 0 --db ${multi_species_db} --report ${sample_id}.${multi_species_db}.multi_profile.tsv ${r1} > ${sample_id}.${multi_species_db}.reads_classification.tsv
    """
    else
    """
    kraken2 --threads ${task.cpus} --confidence 0 --db ${multi_species_db} --report ${sample_id}.${multi_species_db}.multi_profile.tsv --paired ${r1} ${r2} > ${sample_id}.${multi_species_db}.reads_classification.tsv
    """
}

process extract_candidate_species {
    input:
    tuple path(multi_profile), val(taxonomy), val(sample_id), path(r1), path(r2)
    output:
    tuple path("candidates.tsv"), val(taxonomy), val(sample_id), path(r1), path(r2)
    script:
    if (taxonomy == 'gtdb')
    """
    cat $multi_profile | tr -s " " | grep -P "\tS\t" | sort -nr | head -n $params.nspecies | cut -f6 | sed 's/ /_/g' | sed 's/^_//' > candidates.tsv
    """ 
    else
    """
    cat $multi_profile | tr -s " " | grep -P "\tS\t" | sort -nr | head -n $params.nspecies | cut -f5 > candidates.tsv
    """ 
}

process build_species_database {
    label 'medium_task'
    publishDir "species_db/${taxonomy}/", mode: params.publish_mode, overwrite: true
    conda 'conda/kraken.yaml'
    input:
    tuple val(species), val(taxonomy)
    each path(genome_clusters)
    each path(metadata)
    each path(ncbi_dmp)
    each path(gtdb_dmp)
    output:
    path(species)
    script:
    """
    mkdir -p $species/taxonomy
    if [[ $taxonomy == "gtdb" ]]; then
    cp $gtdb_dmp/taxonomy/*.dmp $species/taxonomy
    cat $metadata | cut -f20,110 | sed 's/ /_/g' | grep ${species} | cut -f2 | sort | uniq > groups
    else
    cp $ncbi_dmp/taxonomy/*.dmp $species/taxonomy
    cat $metadata | cut -f80,110 | sed 's/ /_/g' | grep ${species} | cut -f2 | sort | uniq > groups
    fi
    for group in \$(cat groups)
    do
    find $genome_clusters/\$group/ -name "*.selected.${taxonomy}.fna" -exec kraken2-build --add-to-library {} --db ${species} \\; || true
    done
    kraken2-build --build --db ${species} --threads $task.cpus
    kraken2-build --clean --db ${species}
    """
}

process single_profiling {
    label 'medium_task'
    publishDir "results_profiling/single/reads/", mode: params.publish_mode, overwrite: true, pattern: '*.reads_classification.tsv'
    publishDir "results_profiling/single/profile/", mode: params.publish_mode, overwrite: true, pattern: '*.single_profile.tsv'
    conda 'conda/kraken.yaml'
    input:
    tuple val(species), val(taxonomy), val(sample_id), path(r1), path(r2)
    path(all_single_dbs)
    output:
    tuple path("${sample_id}.${taxonomy}.${species}.single_profile.tsv"), path("${sample_id}.${taxonomy}.${species}.reads_classification.tsv"), val(taxonomy), val(sample_id)
    script:
    if (r2.toString() == 'null')
    """
    kraken2 --threads ${task.cpus} --confidence 0 --db ${species} --report ${sample_id}.${taxonomy}.${species}.single_profile.tsv ${r1} > ${sample_id}.${taxonomy}.${species}.reads_classification.tsv
    """
    else
    """
    kraken2 --threads ${task.cpus} --confidence 0 --db ${species} --report ${sample_id}.${taxonomy}.${species}.single_profile.tsv --paired ${r1} ${r2} > ${sample_id}.${taxonomy}.${species}.reads_classification.tsv
    """
}

process taxonomic_confidence {
    label 'medium_task'
    publishDir "confidence_plots", mode: params.publish_mode, overwrite: true
    conda 'conda/confidence.yaml'
    input:
    each path(taxonomic_confidence_script)
    tuple path(profile), path(reads), val(taxonomy), val(sample_id)
    output:
    path("*.confidence_plot.png")
    script:
    """
    python3 ${taxonomic_confidence_script} ${params.nspecies} ${task.cpus} ${sample_id} ${taxonomy}
    """
}

workflow {
    // structure the single samples in a paired-like channel, the null value will be taken into account by the processes running kraken
    paired_samples.concat(single_samples.map(it -> [it.toString().split('/')[-1].toString().split('\\.')[0], it, "/dev/null"])).set{ samples }
    samples.view{ it -> "processing ${it[0]}\n" }

    if (params.ncbi & params.gtdb) {
        ncbi_multi_species_db.concat(gtdb_multi_species_db).set { multi_species_dbs }
    }
    else if (params.ncbi) {
        channel.of("ncbi").set{ taxonomy }
        ncbi_multi_species_db.set { multi_species_dbs }
    }
    else if (params.gtdb) {
        channel.of("gtdb").set{ taxonomy }
        gtdb_multi_species_db.set { multi_species_dbs }
    }
    // run with the multi-species db (ncbi and/or gtdb), add profiles info to the outputss
    multi_profiling(samples, multi_species_dbs).set{ samples_with_multi_profiles_and_reads_classification }

    // identifiy the n species to process with single species db
    // convert the taxonoy info into actual String
    // it[0] is multi profile
    // 1 is taxonomy info
    // 2-4 for keeping track of sample
    // 5 reads classification is dropped
    extract_candidate_species(samples_with_multi_profiles_and_reads_classification.map{it -> [it[0], it[1].toString(), it[2], it[3], it[4]]}).set{ samples_with_candidate_species }

    // build the species single db, need to convert a list of one species into a string after splitCsv (it[0][0])
    // it[0] is candidate species
    // 1 is taxonomy info
    build_species_database(samples_with_candidate_species.splitCsv().map{ it -> [it[0][0], it[1]] }.unique(), genome_clusters, metadata, ncbi_multi_species_db, gtdb_multi_species_db).set{ single_species_dbs}

    // run with the species single dbs
    // it[0] is candidate species
    // 1 is taxonomy info
    // 2-4 for keeping track of sample
    single_species_dbs.unique().collect().set{ all_single_dbs }
    single_profiling(samples_with_candidate_species.splitCsv().map{ it -> [it[0][0], it[1], it[2], it[3], it[4]] }, all_single_dbs).set{ samples_with_single_profiles_and_reads_classification }

    // produce the confidence plots by exploring the read classification with different confidence values
    taxonomic_confidence(taxonomic_confidence_script, samples_with_single_profiles_and_reads_classification.groupTuple(by: [2,3]))

}
