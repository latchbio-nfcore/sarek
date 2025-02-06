from typing import List, Optional

from latch.resources.launch_plan import LaunchPlan
from latch.resources.workflow import workflow
from latch.types import metadata
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile

from wf.entrypoint import (
    Aligner,
    GenomeReference,
    ReferenceType,
    Sample,
    StepOptions,
    initialize,
    nextflow_runtime,
)


@workflow(metadata._nextflow_metadata)
def nf_nf_core_sarek(
    run_name: str,
    input: List[Sample],
    genome_source: str,
    wes: bool,
    intervals: Optional[LatchFile],
    no_intervals: bool,
    tools: Optional[str],
    skip_tools: Optional[str],
    trim_fastq: bool,
    umi_read_structure: Optional[str],
    save_mapped: bool,
    save_output_as_bam: bool,
    use_gatk_spark: Optional[str],
    concatenate_vcfs: bool,
    only_paired_variant_calling: bool,
    joint_germline: bool,
    joint_mutect2: bool,
    bcftools_annotations: Optional[LatchFile],
    bcftools_annotations_tbi: Optional[LatchFile],
    bcftools_header_lines: Optional[LatchFile],
    dbsnp_vqsr: Optional[str],
    fasta: Optional[LatchFile],
    fasta_fai: Optional[LatchFile],
    known_indels_vqsr: Optional[str],
    known_indels: Optional[LatchFile],
    dbsnp: Optional[LatchFile],
    known_snps: Optional[LatchFile],
    known_snps_tbi: Optional[LatchFile],
    known_snps_vqsr: Optional[str],
    ngscheckmate_bed: Optional[LatchFile],
    snpeff_db: Optional[str],
    snpeff_genome: Optional[str],
    vep_genome: Optional[str],
    vep_species: Optional[str],
    vep_cache_version: Optional[str],
    save_reference: bool,
    build_only_index: bool,
    download_cache: bool,
    igenomes_base: Optional[LatchDir],
    igenomes_ignore: bool,
    vep_cache: Optional[LatchDir],
    snpeff_cache: Optional[LatchDir],
    email: Optional[str],
    multiqc_title: Optional[str],
    genome: Optional[GenomeReference],
    multiqc_methods_description: Optional[str],
    # latch_genome: ReferenceType = ReferenceType.homo_sapiens,
    step: StepOptions = StepOptions.mapping,
    split_fastq: int = 50000000,
    nucleotides_per_second: int = 200000,
    aligner: Aligner = Aligner.bwa_mem,
    vep_custom_args: Optional[
        str
    ] = "--everything --filter_common --per_gene --total_length --offline --format vcf",
    vep_version: Optional[str] = "111.0-0",
    outdir: LatchOutputDir = LatchOutputDir("latch:///Sarek"),
) -> None:
    """
    nf-core/sarek is a workflow designed to detect variants on whole genome or targeted sequencing data. It can also handle tumour / normal pairs and could include additional relapses.

    <html>
    <p align="center">
    <img src="https://user-images.githubusercontent.com/31255434/182289305-4cc620e3-86ae-480f-9b61-6ca83283caa5.jpg" alt="Latch Verified" width="100">
    </p>

    <p align="center">
    <strong>
    Latch Verified
    </strong>
    </p>

    [![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.3476425-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.3476425)

    ## Introduction

    **nf-core/sarek** is a workflow designed to detect variants on whole genome or targeted sequencing data. Initially designed for Human, and Mouse, it can work on any species with a reference genome. Sarek can also handle tumour / normal pairs and could include additional relapses.

    The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

    This workflow is hosted on Latch Workflows, using a native Nextflow integration, with a graphical interface for accessible analysis by scientists. There is also an integration with Latch Registry so that batched workflows can be launched from “graphical sample sheets” or tables associating raw sequencing files with metadata.

    ## Pipeline summary

    Depending on the options and samples provided, the pipeline can currently perform the following:

    - Form consensus reads from UMI sequences (`fgbio`)
    - Sequencing quality control and trimming (enabled by `--trim_fastq`) (`FastQC`, `fastp`)
    - Map Reads to Reference (`BWA-mem`, `BWA-mem2`, `dragmap` or `Sentieon BWA-mem`)
    - Process BAM file (`GATK MarkDuplicates`, `GATK BaseRecalibrator` and `GATK ApplyBQSR` or `Sentieon LocusCollector` and `Sentieon Dedup`)
    - Summarise alignment statistics (`samtools stats`, `mosdepth`)
    - Variant calling (enabled by `--tools`, see [compatibility](https://nf-co.re/sarek/latest/docs/usage#which-variant-calling-tool-is-implemented-for-which-data-type)):
    - `ASCAT`
    - `CNVkit`
    - `Control-FREEC`
    - `DeepVariant`
    - `freebayes`
    - `GATK HaplotypeCaller`
    - `Manta`
    - `mpileup`
    - `MSIsensor-pro`
    - `Mutect2`
    - `Sentieon Haplotyper`
    - `Strelka2`
    - `TIDDIT`
    - Variant filtering and annotation (`SnpEff`, `Ensembl VEP`, `BCFtools annotate`)
    - Summarise and represent QC (`MultiQC`)

    ## Credits

    Sarek was originally written by Maxime U Garcia and Szilveszter Juhos at the [National Genomics Infastructure](https://ngisweden.scilifelab.se) and [National Bioinformatics Infastructure Sweden](https://nbis.se) which are both platforms at [SciLifeLab](https://scilifelab.se), with the support of [The Swedish Childhood Tumor Biobank (Barntumörbanken)](https://ki.se/forskning/barntumorbanken).
    Friederike Hanssen and Gisela Gabernet at [QBiC](https://www.qbic.uni-tuebingen.de/) later joined and helped with further development.

    The Nextflow DSL2 conversion of the pipeline was lead by Friederike Hanssen and Maxime U Garcia.

    Maintenance is now lead by Friederike Hanssen and Maxime U Garcia (now at [Seqera Labs](https://seqera/io))

    Main developers:

    - [Maxime U Garcia](https://github.com/maxulysse)
    - [Friederike Hanssen](https://github.com/FriederikeHanssen)

    We thank the following people for their extensive assistance in the development of this pipeline:

    - [Abhinav Sharma](https://github.com/abhi18av)
    - [Adam Talbot](https://github.com/adamrtalbot)
    - [Adrian Lärkeryd](https://github.com/adrlar)
    - [Alexander Peltzer](https://github.com/apeltzer)
    - [Alison Meynert](https://github.com/ameynert)
    - [Anders Sune Pedersen](https://github.com/asp8200)
    - [arontommi](https://github.com/arontommi)
    - [BarryDigby](https://github.com/BarryDigby)
    - [Bekir Ergüner](https://github.com/berguner)
    - [bjornnystedt](https://github.com/bjornnystedt)
    - [cgpu](https://github.com/cgpu)
    - [Chela James](https://github.com/chelauk)
    - [David Mas-Ponte](https://github.com/davidmasp)
    - [Edmund Miller](https://github.com/edmundmiller)
    - [Famke Bäuerle](https://github.com/famosab)
    - [Francesco Lescai](https://github.com/lescai)
    - [Gavin Mackenzie](https://github.com/GCJMackenzie)
    - [Gisela Gabernet](https://github.com/ggabernet)
    - [Grant Neilson](https://github.com/grantn5)
    - [gulfshores](https://github.com/gulfshores)
    - [Harshil Patel](https://github.com/drpatelh)
    - [Hongwei Ye](https://github.com/YeHW)
    - [James A. Fellows Yates](https://github.com/jfy133)
    - [Jesper Eisfeldt](https://github.com/J35P312)
    - [Johannes Alneberg](https://github.com/alneberg)
    - [José Fernández Navarro](https://github.com/jfnavarro)
    - [Júlia Mir Pedrol](https://github.com/mirpedrol)
    - [Ken Brewer](https://github.com/kenibrewer)
    - [Lasse Westergaard Folkersen](https://github.com/lassefolkersen)
    - [Lucia Conde](https://github.com/lconde-ucl)
    - [Malin Larsson](https://github.com/malinlarsson)
    - [Marcel Martin](https://github.com/marcelm)
    - [Nick Smith](https://github.com/nickhsmith)
    - [Nicolas Schcolnicov](https://github.com/nschcolnicov)
    - [Nilesh Tawari](https://github.com/nilesh-tawari)
    - [Nils Homer](https://github.com/nh13)
    - [Olga Botvinnik](https://github.com/olgabot)
    - [Oskar Wacker](https://github.com/WackerO)
    - [pallolason](https://github.com/pallolason)
    - [Paul Cantalupo](https://github.com/pcantalupo)
    - [Phil Ewels](https://github.com/ewels)
    - [Sabrina Krakau](https://github.com/skrakau)
    - [Sam Minot](https://github.com/sminot)
    - [Sebastian-D](https://github.com/Sebastian-D)
    - [Silvia Morini](https://github.com/silviamorins)
    - [Simon Pearce](https://github.com/SPPearce)
    - [Solenne Correard](https://github.com/scorreard)
    - [Susanne Jodoin](https://github.com/SusiJo)
    - [Szilveszter Juhos](https://github.com/szilvajuhos)
    - [Tobias Koch](https://github.com/KochTobi)
    - [Winni Kretzschmar](https://github.com/winni2k)

    ## Acknowledgements

    |      [![Barntumörbanken](docs/images/BTB_logo.png)](https://ki.se/forskning/barntumorbanken)      |            [![SciLifeLab](docs/images/SciLifeLab_logo.png)](https://scilifelab.se)             |
    | :-----------------------------------------------------------------------------------------------: | :--------------------------------------------------------------------------------------------: |
    | [![National Genomics Infrastructure](docs/images/NGI_logo.png)](https://ngisweden.scilifelab.se/) | [![National Bioinformatics Infrastructure Sweden](docs/images/NBIS_logo.png)](https://nbis.se) |
    |              [![QBiC](docs/images/QBiC_logo.png)](https://www.qbic.uni-tuebingen.de)              |                   [![GHGA](docs/images/GHGA_logo.png)](https://www.ghga.de/)                   |
    |                     [![DNGC](docs/images/DNGC_logo.png)](https://eng.ngc.dk/)                     |                                                                                                |

    ## Contributions & Support

    If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

    For further information or help, don't hesitate to get in touch on the [Slack `#sarek` channel](https://nfcore.slack.com/channels/sarek) (you can join with [this invite](https://nf-co.re/join/slack)), or contact us: [Maxime U Garcia](mailto:maxime.garcia@seqera.io?subject=[GitHub]%20nf-core/sarek), [Friederike Hanssen](mailto:friederike.hanssen@qbic.uni-tuebingen.de?subject=[GitHub]%20nf-core/sarek)

    ## Citations

    If you use `nf-core/sarek` for your analysis, please cite the `Sarek` article as follows:

    > Friederike Hanssen, Maxime U Garcia, Lasse Folkersen, Anders Sune Pedersen, Francesco Lescai, Susanne Jodoin, Edmund Miller, Oskar Wacker, Nicholas Smith, nf-core community, Gisela Gabernet, Sven Nahnsen **Scalable and efficient DNA sequencing analysis on different compute infrastructures aiding variant discovery** _NAR Genomics and Bioinformatics_ Volume 6, Issue 2, June 2024, lqae031, [doi: 10.1093/nargab/lqae031](https://doi.org/10.1093/nargab/lqae031).

    > Garcia M, Juhos S, Larsson M et al. **Sarek: A portable workflow for whole-genome sequencing analysis of germline and somatic variants [version 2; peer review: 2 approved]** _F1000Research_ 2020, 9:63 [doi: 10.12688/f1000research.16665.2](http://dx.doi.org/10.12688/f1000research.16665.2).

    You can cite the sarek zenodo record for a specific version using the following [doi: 10.5281/zenodo.3476425](https://doi.org/10.5281/zenodo.3476425)

    An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

    You can cite the `nf-core` publication as follows:

    > **The nf-core framework for community-curated bioinformatics pipelines.**
    >
    > Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
    >
    > _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

    """

    pvc_name: str = initialize(run_name=run_name)
    nextflow_runtime(
        run_name=run_name,
        pvc_name=pvc_name,
        input=input,
        step=step,
        genome_source=genome_source,
        # latch_genome=latch_genome,
        outdir=outdir,
        split_fastq=split_fastq,
        wes=wes,
        intervals=intervals,
        nucleotides_per_second=nucleotides_per_second,
        no_intervals=no_intervals,
        tools=tools,
        skip_tools=skip_tools,
        trim_fastq=trim_fastq,
        umi_read_structure=umi_read_structure,
        aligner=aligner,
        save_mapped=save_mapped,
        save_output_as_bam=save_output_as_bam,
        use_gatk_spark=use_gatk_spark,
        concatenate_vcfs=concatenate_vcfs,
        only_paired_variant_calling=only_paired_variant_calling,
        joint_germline=joint_germline,
        joint_mutect2=joint_mutect2,
        vep_custom_args=vep_custom_args,
        vep_version=vep_version,
        bcftools_annotations=bcftools_annotations,
        bcftools_annotations_tbi=bcftools_annotations_tbi,
        bcftools_header_lines=bcftools_header_lines,
        genome=genome,
        dbsnp_vqsr=dbsnp_vqsr,
        fasta=fasta,
        fasta_fai=fasta_fai,
        known_indels_vqsr=known_indels_vqsr,
        known_indels=known_indels,
        dbsnp=dbsnp,
        known_snps=known_snps,
        known_snps_tbi=known_snps_tbi,
        known_snps_vqsr=known_snps_vqsr,
        ngscheckmate_bed=ngscheckmate_bed,
        snpeff_db=snpeff_db,
        snpeff_genome=snpeff_genome,
        vep_genome=vep_genome,
        vep_species=vep_species,
        vep_cache_version=vep_cache_version,
        save_reference=save_reference,
        build_only_index=build_only_index,
        download_cache=download_cache,
        igenomes_base=igenomes_base,
        igenomes_ignore=igenomes_ignore,
        vep_cache=vep_cache,
        snpeff_cache=snpeff_cache,
        email=email,
        multiqc_title=multiqc_title,
        multiqc_methods_description=multiqc_methods_description,
    )


LaunchPlan(
    nf_nf_core_sarek,
    "Test Data",
    {
        "input": [
            Sample(
                patient="TestSample",
                sex="XX",
                status=0,
                sample="TestSample",
                lane=1,
                fastq_1=LatchFile(
                    "s3://latch-public/nf-core/sarek/test_data/test_L1_1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/nf-core/sarek/test_data/test_L1_2.fastq.gz"
                ),
                bam=None,
                bai=None,
                cram=None,
                crai=None,
                table=None,
                vcf=None,
            ),
            Sample(
                patient="TestSample",
                sex="XX",
                status=0,
                sample="TestSample",
                lane=2,
                fastq_1=LatchFile(
                    "s3://latch-public/nf-core/sarek/test_data/test_L2_1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/nf-core/sarek/test_data/test_L2_2.fastq.gz"
                ),
                bam=None,
                bai=None,
                cram=None,
                crai=None,
                table=None,
                vcf=None,
            ),
        ],
        "run_name": "Test_Run",
        "tools": "strelka",
        # "trim_fastq": False,
        "bcftools_annotations": LatchFile(
            "s3://latch-public/nf-core/sarek/test_data/test2.vcf.gz"
        ),
        "bcftools_annotations_tbi": LatchFile(
            "s3://latch-public/nf-core/sarek/test_data/test2.vcf.gz.tbi"
        ),
        "bcftools_header_lines": LatchFile(
            "s3://latch-public/nf-core/sarek/test_data/bcfann_test_header.txt"
        ),
        "genome": GenomeReference.GATK_GRCh37,
        # "step": StepOptions.mapping,
        "split_fastq": 0,
        # "aligner": Aligner.bwa_mem,
    },
)
