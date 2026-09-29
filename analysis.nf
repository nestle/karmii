/*
 * See README.md for instructions on how to run the pipeline and the expected outputs.
 * See LICENSE.md and CONTRIBUTING.md for license and contribution details.
 */

 /*
 * This nf pipeline performs taxonomic profiling of metagenomic samples using Kraken2.
 * It takes as input the raw reads of the samples,
 * the Kraken2 databases, 
 * both for NCBI and GTDB taxonomies,
 * the genome clusters for building the single-species databases.
 * It performs a first round of profiling using the multi-species databases,
 * extracts candidate species from the profiles,
 * and builds single-species databases for the identified candidates.
 * It then performs a second round of profiling using the single-species databases,
 * and generates confidence plots for the taxonomic assignments. 
 */

/*
 * Processes
 */
process multi_profiling {
    label 'big_task'
    publishDir { "results_profiling/multi/reads/${sample_id}/" },
        mode: params.publish_mode,
        overwrite: true,
        pattern: '*.reads_classification.tsv'
    publishDir { "results_profiling/multi/profile/${sample_id}/" },
        mode: params.publish_mode,
        overwrite: true,
        pattern: '*.multi_profile.tsv'
    conda 'conda/kraken.yaml'

    input:
    tuple val(sample_id), path(r1), path(r2)
    each path(multi_species_db)

    output:
    tuple path("${sample_id}.${multi_species_db}.multi_profile.tsv"),
        val(multi_species_db),
        val(sample_id),
        path(r1),
        path(r2),
        path("${sample_id}.${multi_species_db}.reads_classification.tsv")

    script:
    if (r2.toString() == 'null') {
        """
        kraken2 --threads ${task.cpus} --confidence 0 --db ${multi_species_db} \
            --report ${sample_id}.${multi_species_db}.multi_profile.tsv ${r1} \
            > ${sample_id}.${multi_species_db}.reads_classification.tsv
        """
    } else {
        """
        kraken2 --threads ${task.cpus} --confidence 0 --db ${multi_species_db} \
            --report ${sample_id}.${multi_species_db}.multi_profile.tsv \
            --paired ${r1} ${r2} \
            > ${sample_id}.${multi_species_db}.reads_classification.tsv
        """
    }
}

process extract_candidate_species {
    input:
    tuple path(multi_profile), val(taxonomy), val(sample_id), path(r1), path(r2)

    output:
    tuple path('candidates.tsv'), val(taxonomy), val(sample_id), path(r1), path(r2)

    script:
    if (taxonomy == 'gtdb') {
        """
        for genus in \$(cat $multi_profile | tr -s " " | grep -P "\tG\t" | sort -nr | head -n $params.n_genus_species | cut -f6);do \
        cat $multi_profile | tr -s " " | sort -nr | grep -P "\tS\t" | grep -P " \$(echo \$genus | sed 's/g_/s_/g') " | \
        head -n $params.n_genus_species;done | cut -f6 | tr -s " " | sed 's/ /_/g' | sed 's/^_//' > candidates.tsv
        """
    } else {
        """
        for genus in \$(cat $multi_profile | tr -s " " | grep -P "\tG\t" | sort -nr | head -n $params.n_genus_species | cut -f6);do \
        cat $multi_profile | tr -s " " | sort -nr | grep -P "\tS\t" | grep -P " \$genus " | \
        head -n $params.n_genus_species;done | cut -f5 > candidates.tsv
        """
    }
}

process build_species_database {
    label 'medium_task'
    publishDir   { "species_db/${taxonomy}/" }, mode: params.publish_mode, overwrite: true
    conda 'conda/kraken.yaml'

    input:
    tuple val(species), val(taxonomy)
    each path(genome_clusters)
    each path(metadata)
    each path(ncbi_multi_species_db)
    each path(gtdb_multi_species_db)

    output:
    path(species)

    script:
    """
    mkdir -p $species/taxonomy
    if [[ $taxonomy == "gtdb" ]]; then
        cp $gtdb_multi_species_db/taxonomy/*.dmp $species/taxonomy
        cat $metadata | cut -f20,110 | sed 's/ /_/g' | grep -P "${species}\t" | cut -f2 | \
            sort | uniq > groups
    else
        cp $ncbi_multi_species_db/taxonomy/*.dmp $species/taxonomy
        cat $metadata | cut -f76,110 | sed 's/ /_/g' | grep -P "^${species}\t" | cut -f2 | \
            sort | uniq > groups
    fi

    for group in \$(cat groups)
    do
        find $genome_clusters/\$group/ -name "*.selected.${taxonomy}.fna" \
            -exec kraken2-build --no-masking --add-to-library {} --db ${species} \\; || true
    done

    kraken2-build --build --db ${species} --threads $task.cpus
    kraken2-build --clean --db ${species}
    """
}

process single_profiling {
    label 'medium_task'
    publishDir { "results_profiling/single/reads/${sample_id}/" },
        mode: params.publish_mode,
        overwrite: true,
        pattern: '*.reads_classification.tsv'
    publishDir { "results_profiling/single/profile/${sample_id}/" },
        mode: params.publish_mode,
        overwrite: true,
        pattern: '*.single_profile.tsv'
    conda 'conda/kraken.yaml'

    input:
    tuple val(species), val(taxonomy), val(sample_id), path(r1), path(r2), path(db)

    output:
    tuple path("${sample_id}.${taxonomy}.${species}.single_profile.tsv"),
        path("${sample_id}.${taxonomy}.${species}.reads_classification.tsv"),
        val(taxonomy),
        val(sample_id)

    script:
    if (r2.toString() == 'null') {
        """
        kraken2 --threads ${task.cpus} --confidence 0 --db ${species} \
            --report ${sample_id}.${taxonomy}.${species}.single_profile.tsv ${r1} \
            > ${sample_id}.${taxonomy}.${species}.reads_classification.tsv
        """
    } else {
        """
        kraken2 --threads ${task.cpus} --confidence 0 --db ${species} \
            --report ${sample_id}.${taxonomy}.${species}.single_profile.tsv \
            --paired ${r1} ${r2} \
            > ${sample_id}.${taxonomy}.${species}.reads_classification.tsv
        """
    }
}

process taxonomic_confidence {
    label 'medium_task'
    publishDir { "confidence_plots/${sample_id}/" }, mode: params.publish_mode, overwrite: true, pattern: '*.confidence_plot.*'
    conda 'conda/confidence.yaml'

    input:
    each path(taxonomic_confidence_script)
    tuple path(profile), path(reads), val(taxonomy), val(sample_id)

    output:
    tuple path('*.confidence_plot.png'), path('*.confidence_plot.svg'), path('values_for_decision.txt'), val(taxonomy), val(sample_id)

    script:
    """
    python3 ${taxonomic_confidence_script} \$(ls *.single_profile.tsv | wc -l) ${task.cpus} ${sample_id} ${taxonomy}
    # and extract the important values for the summary step
    for x in *.*.*.reads_classification.tsv.conf.txt;do a=\$(grep -P "^0.5000," \$x | cut -f2 -d",");b=\$(grep -P "^0.0500," \$x | cut -f2 -d",");echo \$a","\$b",";done | sort -nr | head -n 2 > values_for_decision.txt
    """
}

process summary {
    label 'medium_task'
    publishDir { "./tmp" }, mode: params.publish_mode, overwrite: true
    input:
    tuple val(taxonomy), val(sample_id), path(single_values_for_decision), path(multi_profiles)

    output:
    path('*.summary.tsv')

    script:
    """
    #!/usr/bin/env python
    r=[[float(x) for x in l.strip().split(',')[0:2]] for l in open('${single_values_for_decision}')]
    # r[0] is first species, r[1] is second species
    # r[x][0] is the confidence at 0.5
    # r[x][1] is the confidence at 0.05
    with open('${sample_id}.${taxonomy}.summary.tsv', 'w') as out:
        try:
            # validate that the slope is above -20%
            test1=(r[0][0]-r[0][1])/(0.5-0.05) > -0.2
            # validate that the confidence at 0.05 is above 0.95
            test2=r[0][1]>0.95
            # validate the separation of the two lines at conf 0.5
            test3=r[1][0]<=r[0][0]-0.05
            out.write(f'${sample_id}\t${taxonomy}\t{test1}\t{test2}\t{test3}\t{"PASS" if test1 and test2 and test3 else "FAIL"}\\n')
        except IndexError:
            out.write(f'${sample_id}\t${taxonomy}\tERROR\tERROR\tERROR\tFAIL\\n')
        """
}

process concat_summary {
    label 'medium_task'
    publishDir { "./" }, mode: params.publish_mode, overwrite: true
    input:
    path(summary)

    output:
    path('summary.tsv')

    script:
    """
    echo -e 'sample_id\ttaxonomy\tslope > -20%\ty@x=0.05 > 0.95\tsp2@0.5 << sp1@0.5\tstatus' > summary.tsv
    cat ${summary} | sort -n >> summary.tsv
    """
}

/*
 * Workflow
 */
workflow {

    /*
    * channels
    */
    channel
        .fromFilePairs('samples/*_R{1,2}*.fastq.gz', flat: true)
        .set { paired_samples_gz }

    channel
        .fromFilePairs('samples/*_R{1,2}*.fastq', flat: true)
        .set { paired_samples }

    channel
        .fromPath('samples/*.fastq.gz')
        .filter { file ->
            !(file.name =~ /.*_R1.*\.fastq.gz$/ || file.name =~ /.*_R2.*\.fastq.gz$/)
        }
        .set { single_samples_gz }

    channel
        .fromPath('samples/*.fastq')
        .filter { file ->
            !(file.name =~ /.*_R1.*\.fastq$/ || file.name =~ /.*_R2.*\.fastq$/)
        }
        .set { single_samples }

    channel
        .fromPath('samples/*.fasta.gz')
        .set { fasta_samples_gz }

    channel
        .fromPath('samples/*.fasta')
        .set { fasta_samples }

    channel
        .fromPath('krakendbs/ncbi/')
        .set { ncbi_multi_species_db }

    channel
        .fromPath('krakendbs/gtdb/')
        .set { gtdb_multi_species_db }

    channel
        .fromPath('clusters/')
        .set { genome_clusters }

    channel
        .fromPath('metadata/metadata_final.tsv')
        .set { metadata }

    channel
        .fromPath('scripts/taxonomic_confidence.py')
        .set { taxonomic_confidence_script }

    // fuse the samples channels and
    // add a placeholder for the second read
    // in the case of single-end samples
    paired_samples_gz
        .concat(paired_samples)
        .concat(
            single_samples_gz.map { it ->
                [
                    it.toString().split('/')[-1].toString().split('\\.')[0],
                    it,
                    '/dev/null'
                ]
            }
        )
        .concat(
            single_samples.map { it ->
                [
                    it.toString().split('/')[-1].toString().split('\\.')[0],
                    it,
                    '/dev/null'
                ]
            }
        )
        .concat(
            fasta_samples_gz.map { it ->
                [
                    it.toString().split('/')[-1].toString().split('\\.')[0],
                    it,
                    '/dev/null'
                ]
            }
        )
        .concat(
            fasta_samples.map { it ->
                [
                    it.toString().split('/')[-1].toString().split('\\.')[0],
                    it,
                    '/dev/null'
                ]
            }
        )
        .set { samples }

    if (params.ncbi & params.gtdb) {
        ncbi_multi_species_db
            .concat(gtdb_multi_species_db)
            .set { multi_species_dbs }
    } else if (params.ncbi) {
        channel.of('ncbi').set { taxonomy }
        ncbi_multi_species_db.set { multi_species_dbs }
    } else if (params.gtdb) {
        channel.of('gtdb').set { taxonomy }
        gtdb_multi_species_db.set { multi_species_dbs }
    }

    multi_profiling(samples, multi_species_dbs)
        .set { samples_with_multi_profiles_and_reads_classification }

    extract_candidate_species(
        samples_with_multi_profiles_and_reads_classification.map { it ->
            [it[0], it[1].toString(), it[2], it[3], it[4]]
        }
    ).set { samples_with_candidate_species }

    build_species_database(
        samples_with_candidate_species
            .splitCsv()
            .map { it -> [it[0][0], it[1]] }
            .unique(),
        genome_clusters,
        metadata,
        ncbi_multi_species_db,
        gtdb_multi_species_db
    ).set { single_species_dbs }

    single_species_dbs
        .map { p -> [p.name.split('/')[-1], p] }
        .set { dbs_with_species_as_key } // to be used below for combining with the samples channel

    samples_with_candidate_species
        .splitCsv()
        .map { it -> [it[0][0], it[1], it[2], it[3], it[4]] }  // [species, taxonomy, sample_id, r1, r2]
        .combine(dbs_with_species_as_key, by: 0)             // combine to provide the corresponding database → [species, taxonomy, sample_id, r1, r2, db_path]
        .set { samples_for_single_profiling }
    
    single_profiling(
        samples_for_single_profiling
    ).set { samples_with_single_profiles_and_reads_classification }

    taxonomic_confidence(
        taxonomic_confidence_script,
        samples_with_single_profiles_and_reads_classification.groupTuple(by: [2, 3]) // group by taxonomy and sample_id → [profile, reads, taxonomy, sample_id]
    ).set { taxonomic_confidence_results }

    // produce a summary based on the multi profile and the single profiles important values
    samples_with_multi_profiles_and_reads_classification.map{it -> [it[0],it[1],it[2]]}.set{ multi  }
    taxonomic_confidence_results.map{ it -> [it[2], it[3], it[4]] }.set{ single_values_for_decision }
    concat_summary(summary(single_values_for_decision.join(multi,by:[2,3]).map{ it -> [it[2], it[0], it[1], it[3]] }).collect())

}
