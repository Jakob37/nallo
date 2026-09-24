process GENERATE_GENS_METADATA {
    tag "${meta.id}"
    label 'process_single'

    conda 'conda-forge::python=3.8.3'
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.8.3'
        : 'biocontainers/python:3.8.3'}"

    input:
    tuple val(meta), path(ratios), path(roh), path(upd_regions), path(upd_sites)
    path fai
    path metadata_script

    output:
    tuple val(meta), path("${meta.id}.meta.tsv"), emit: sample_meta
    tuple val(meta), path("${meta.id}.chrom_meta.tsv"), emit: chrom_meta
    tuple val(meta), path("${meta.id}.gens_track.roh.bed"), emit: roh_track
    tuple val(meta), path("${meta.id}.gens_track.upd.bed"), optional: true, emit: upd_track
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    script:
    if (!(meta.sex in [1, 2])) {
        error("Gens metadata requires sex 1 or 2 for ${meta.id}")
    }
    def sex = meta.sex == 1 ? 'M' : 'F'
    def upd_args = upd_regions && upd_sites ? "--upd-regions ${upd_regions} --upd-sites ${upd_sites}" : ''
    """
    python ${metadata_script} \\
        --sample '${meta.id}' \\
        --sex ${sex} \\
        --fai ${fai} \\
        --ratios ${ratios} \\
        --roh ${roh} \\
        ${upd_args}
    """

    stub:
    """
    touch ${meta.id}.meta.tsv ${meta.id}.chrom_meta.tsv ${meta.id}.gens_track.roh.bed
    ${upd_regions && upd_sites ? "touch ${meta.id}.gens_track.upd.bed" : 'true'}
    """
}
