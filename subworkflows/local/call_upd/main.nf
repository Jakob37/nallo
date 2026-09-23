include { BCFTOOLS_VIEW as CONVERT_UPD_BCF            } from '../../../modules/nf-core/bcftools/view/main'
include { UPD as UPD_REGIONS                          } from '../../../modules/nf-core/upd/main'
include { UPD as UPD_SITES                            } from '../../../modules/nf-core/upd/main'
include { VCF_CONCAT_SORT_VARIANTS as CONCAT_SORT_UPD } from '../vcf_concat_sort_variants/main'

workflow CALL_UPD {
    take:
    ch_family_bcf
    ch_samplesheet

    main:
    // Select every child whose named mother and father are present in the family.
    ch_samplesheet
        .map { meta, _reads -> [meta.family_id, meta] }
        .groupTuple()
        .flatMap { family_id, samples -> eligibleUpdTrios(family_id, samples) }
        .set { ch_upd_trios }

    // Convert the unfiltered Echtvar output to indexed VCF, but only for families
    // with a complete trio.
    ch_family_bcf
        .map { meta, bcf -> [meta.family_id, meta, bcf] }
        .combine(ch_upd_trios.map { family_id, _child_meta, _mother, _father -> [family_id] }.unique(), by: 0)
        .map { _family_id, meta, bcf -> [meta, bcf, []] }
        .set { ch_upd_family_bcf }

    CONVERT_UPD_BCF(
        ch_upd_family_bcf,
        [],
        [],
        [],
    )

    // Combine all annotated call regions into one family VCF.
    CONVERT_UPD_BCF.out.vcf
        .join(CONVERT_UPD_BCF.out.tbi, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, vcf, index ->
            [groupKey([id: meta.family_id], meta.num_intervals), vcf, index]
        }
        .groupTuple()
        .map { key, vcfs, indexes -> [key.getGroupTarget(), vcfs, indexes] }
        .set { ch_upd_family_vcfs }

    CONCAT_SORT_UPD(ch_upd_family_vcfs)

    CONCAT_SORT_UPD.out.vcf
        .map { meta, vcf -> [meta.id, vcf] }
        .combine(ch_upd_trios, by: 0)
        .map { _family_id, vcf, child_meta, mother, father ->
            [child_meta + [mother: mother, father: father], vcf]
        }
        .set { ch_upd_input }

    UPD_REGIONS(ch_upd_input)
    UPD_SITES(ch_upd_input)

    emit:
    regions = UPD_REGIONS.out.bed
    sites   = UPD_SITES.out.bed
}

/** Select every child with both named parents present in the same family. */
def eligibleUpdTrios(family_id, samples) {
    def ids = samples.collect { sample -> sample.id } as Set
    samples
        .findAll { child ->
            def mother = child.maternal_id?.toString()
            def father = child.paternal_id?.toString()
            mother && father && mother != '0' && father != '0' && mother != father && mother != child.id && father != child.id && mother in ids && father in ids
        }
        .collect { child ->
            [family_id, [id: child.id, family_id: family_id], child.maternal_id.toString(), child.paternal_id.toString()]
        }
}
