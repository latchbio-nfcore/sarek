//
// PAIRED VARIANT CALLING
//

include { BAM_VARIANT_CALLING_CNVKIT                    } from '../bam_variant_calling_cnvkit/main'
include { BAM_VARIANT_CALLING_FREEBAYES                 } from '../bam_variant_calling_freebayes/main'
include { BAM_VARIANT_CALLING_MPILEUP as MPILEUP_NORMAL } from '../bam_variant_calling_mpileup/main'
include { BAM_VARIANT_CALLING_MPILEUP as MPILEUP_TUMOR  } from '../bam_variant_calling_mpileup/main'
include { BAM_VARIANT_CALLING_SOMATIC_ASCAT             } from '../bam_variant_calling_somatic_ascat/main'
include { BAM_VARIANT_CALLING_SOMATIC_CONTROLFREEC      } from '../bam_variant_calling_somatic_controlfreec/main'
include { BAM_VARIANT_CALLING_SOMATIC_MANTA             } from '../bam_variant_calling_somatic_manta/main'
include { BAM_VARIANT_CALLING_SOMATIC_MUTECT2           } from '../bam_variant_calling_somatic_mutect2/main'
include { BAM_VARIANT_CALLING_SOMATIC_STRELKA           } from '../bam_variant_calling_somatic_strelka/main'
include { BAM_VARIANT_CALLING_SOMATIC_TIDDIT            } from '../bam_variant_calling_somatic_tiddit/main'
include { MSISENSORPRO_MSI_SOMATIC                      } from '../../../modules/nf-core/msisensorpro/msi_somatic/main'

workflow BAM_VARIANT_CALLING_SOMATIC_ALL {
    take:
    tools                         // Mandatory, list of tools to apply
    cram                          // channel: [mandatory] cram
    bwa                           // channel: [optional] bwa
    cf_chrom_len                  // channel: [optional] controlfreec length file
    chr_files
    dbsnp                         // channel: [mandatory] dbsnp
    dbsnp_tbi                     // channel: [mandatory] dbsnp_tbi
    dict                          // channel: [mandatory] dict
    fasta                         // channel: [mandatory] fasta
    fasta_fai                     // channel: [mandatory] fasta_fai
    germline_resource             // channel: [optional]  germline_resource
    germline_resource_tbi         // channel: [optional]  germline_resource_tbi
    intervals                     // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals
    intervals_bed_gz_tbi          // channel: [mandatory] intervals/target regions index zipped and indexed
    intervals_bed_combined        // channel: [mandatory] intervals/target regions in one file unzipped
    mappability
    msisensorpro_scan             // channel: [optional]  msisensorpro_scan
    panel_of_normals              // channel: [optional]  panel_of_normals
    panel_of_normals_tbi          // channel: [optional]  panel_of_normals_tbi
    allele_files                  // channel: [optional]  ascat allele files
    loci_files                    // channel: [optional]  ascat loci files
    gc_file                       // channel: [optional]  ascat gc content file
    rt_file                       // channel: [optional]  ascat rt file

    main:
    versions          = Channel.empty()

    //TODO: Temporary until the if's can be removed and printing to terminal is prevented with "when" in the modules.config
    vcf_freebayes       = Channel.empty()
    vcf_manta           = Channel.empty()
    vcf_strelka         = Channel.empty()
    out_msisensorpro    = Channel.empty()
    vcf_mutect2         = Channel.empty()
    vcf_tiddit          = Channel.empty()

    if (tools.split(',').contains('ascat')) {
        BAM_VARIANT_CALLING_SOMATIC_ASCAT(
            cram,
            allele_files,
            loci_files,
            intervals_bed_combined,
            fasta,
            gc_file,
            rt_file
        )

        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_ASCAT.out.versions)
    }

    // CONTROLFREEC
    if (tools.split(',').contains('controlfreec')) {
        // Remap channels to match module/subworkflow
        cram_normal = cram.map { meta, normal_cram, normal_crai, tumor_cram, tumor_crai -> [ meta, normal_cram, normal_crai ] }
        cram_tumor = cram.map { meta, normal_cram, normal_crai, tumor_cram, tumor_crai -> [ meta, tumor_cram, tumor_crai ] }

        MPILEUP_NORMAL(
            cram_normal,
            // Remap channel to match module/subworkflow
            dict.map{ it -> [ [ id:'dict' ], it ] },
            fasta,
            intervals
        )

        MPILEUP_TUMOR(
            cram_tumor,
            // Remap channel to match module/subworkflow
            dict.map{ it -> [ [ id:'dict' ], it ] },
            fasta,
            intervals
        )

        mpileup_normal = MPILEUP_NORMAL.out.mpileup
        mpileup_tumor = MPILEUP_TUMOR.out.mpileup
            // Remap channel to match module/subworkflow
        mpileup_pair = mpileup_normal.cross(mpileup_tumor).map{ normal, tumor -> [ normal[0], normal[1], tumor[1], [], [], [], [] ] }

        length_file = cf_chrom_len ?: fasta_fai

        BAM_VARIANT_CALLING_SOMATIC_CONTROLFREEC(
            mpileup_pair,
            fasta,
            length_file,
            dbsnp,
            dbsnp_tbi,
            chr_files,
            mappability,
            intervals_bed_combined
        )

        versions = versions.mix(MPILEUP_NORMAL.out.versions)
        versions = versions.mix(MPILEUP_TUMOR.out.versions)
        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_CONTROLFREEC.out.versions)
    }

    // CNVKIT
    if (tools.split(',').contains('cnvkit')) {
        BAM_VARIANT_CALLING_CNVKIT(
            // Remap channel to match module/subworkflow
            cram.map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai -> [ meta, tumor_cram, normal_cram ] },
            fasta,
            fasta_fai,
            intervals_bed_combined,
            []
        )

        versions = versions.mix(BAM_VARIANT_CALLING_CNVKIT.out.versions)
    }

    // FREEBAYES
    if (tools.split(',').contains('freebayes')) {
        BAM_VARIANT_CALLING_FREEBAYES(
            cram,
            // Remap channel to match module/subworkflow
            dict.map{ it -> [ [ id:'dict' ], it ] },
            fasta,
            fasta_fai,
            intervals
        )

        vcf_freebayes = BAM_VARIANT_CALLING_FREEBAYES.out.vcf
        versions   = versions.mix(BAM_VARIANT_CALLING_FREEBAYES.out.versions)
    }

    // MANTA
    if (tools.split(',').contains('manta')) {
        BAM_VARIANT_CALLING_SOMATIC_MANTA(
            cram,
            // Remap channel to match module/subworkflow
            dict.map{ it -> [ [ id:'dict' ], it ] },
            fasta,
            fasta_fai,
            intervals_bed_gz_tbi
        )

        vcf_manta = BAM_VARIANT_CALLING_SOMATIC_MANTA.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_MANTA.out.versions)
    }

    // STRELKA
    if (tools.split(',').contains('strelka')) {
        // Remap channel to match module/subworkflow
        cram_strelka = (tools.split(',').contains('manta')) ?
            cram.join(BAM_VARIANT_CALLING_SOMATIC_MANTA.out.candidate_small_indels_vcf).join(BAM_VARIANT_CALLING_SOMATIC_MANTA.out.candidate_small_indels_vcf_tbi) :
            cram.map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai -> [ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, [], [] ] }

        BAM_VARIANT_CALLING_SOMATIC_STRELKA(
            cram_strelka,
            // Remap channel to match module/subworkflow
            dict.map{ it -> [ [ id:'dict' ], it ] },
            fasta,
            fasta_fai,
            intervals_bed_gz_tbi
        )

        vcf_strelka = Channel.empty().mix(BAM_VARIANT_CALLING_SOMATIC_STRELKA.out.vcf)
        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_STRELKA.out.versions)
    }

    // MSISENSOR
    if (tools.split(',').contains('msisensorpro')) {
        MSISENSORPRO_MSI_SOMATIC(cram.combine(intervals_bed_combined), fasta, msisensorpro_scan)

        versions = versions.mix(MSISENSORPRO_MSI_SOMATIC.out.versions)
        out_msisensorpro = out_msisensorpro.mix(MSISENSORPRO_MSI_SOMATIC.out.output_report)
    }

    // MUTECT2
    if (tools.split(',').contains('mutect2')) {
        BAM_VARIANT_CALLING_SOMATIC_MUTECT2(
            // Remap channel to match module/subworkflow
            cram.map { meta, normal_cram, normal_crai, tumor_cram, tumor_crai -> [ meta, [ normal_cram, tumor_cram ], [ normal_crai, tumor_crai ] ] },
            fasta,
            fasta_fai,
            dict,
            germline_resource,
            germline_resource_tbi,
            panel_of_normals,
            panel_of_normals_tbi,
            intervals
        )

        vcf_mutect2 = BAM_VARIANT_CALLING_SOMATIC_MUTECT2.out.vcf_filtered
        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_MUTECT2.out.versions)
    }

    // TIDDIT
    if (tools.split(',').contains('tiddit')) {
        BAM_VARIANT_CALLING_SOMATIC_TIDDIT(
            // Remap channel to match module/subworkflow
            cram.map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai -> [ meta, normal_cram, normal_crai ] },
            // Remap channel to match module/subworkflow
            cram.map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai -> [ meta, tumor_cram, tumor_crai ] },
            // Remap channel to match module/subworkflow
            fasta.map{ it -> [ [ id:'fasta' ], it ] },
            bwa)
        vcf_tiddit = BAM_VARIANT_CALLING_SOMATIC_TIDDIT.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_TIDDIT.out.versions)
    }

    vcf_all = Channel.empty().mix(
        vcf_freebayes,
        vcf_manta,
        vcf_mutect2,
        vcf_strelka,
        vcf_tiddit
    )

    emit:
    out_msisensorpro
    vcf_all
    vcf_freebayes
    vcf_manta
    vcf_mutect2
    vcf_strelka
    vcf_tiddit

    versions
}
