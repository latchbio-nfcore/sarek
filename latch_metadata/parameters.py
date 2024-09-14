from dataclasses import dataclass
from enum import Enum
from typing import List, Optional

from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile
from latch.types.metadata import (
    Fork,
    ForkBranch,
    LatchRule,
    NextflowParameter,
    Params,
    Section,
    Spoiler,
    Text,
)


@dataclass(frozen=True)
class Sample:
    patient: str
    sex: Optional[str]
    status: Optional[int]
    sample: str
    lane: Optional[int]
    fastq_1: Optional[LatchFile]
    fastq_2: Optional[LatchFile]
    bam: Optional[LatchFile]
    bai: Optional[LatchFile]
    cram: Optional[LatchFile]
    crai: Optional[LatchFile]
    table: Optional[LatchFile]
    vcf: Optional[LatchFile]


class GenomeReference(Enum):
    GATK_GRCh37 = "GATK Human Genome GRCh37 (hg19) - Broad Institute"
    GATK_GRCh38 = "GATK Human Genome GRCh38 (hg38) - Broad Institute"
    Ensembl_GRCh37 = "Ensembl Human Genome GRCh37 (hg19)"
    NCBI_GRCh38 = "NCBI Human Genome GRCh38 (hg38)"
    CHM13 = "UCSC Human Genome CHM13 (T2T Complete)"
    GRCm38 = "Ensembl Mouse Genome GRCm38 (mm10)"
    TAIR10 = "Ensembl Arabidopsis thaliana Genome TAIR10"
    EB2 = "Ensembl Bacillus subtilis Genome EB2"
    UMD3_1 = "Ensembl Cow Genome UMD3.1"
    WBcel235 = "Ensembl Worm (C. elegans) Genome WBcel235"
    CanFam3_1 = "Ensembl Dog Genome CanFam3.1"
    GRCz10 = "Ensembl Zebrafish Genome GRCz10"
    BDGP6 = "Ensembl Drosophila Genome BDGP6"
    EquCab2 = "Ensembl Horse Genome EquCab2"
    EB1 = "Ensembl Escherichia coli Genome EB1"
    Galgal4 = "Ensembl Chicken Genome Galgal4"
    Gm01 = "Ensembl Soybean Genome Gm01"
    Mmul_1 = "Ensembl Rhesus Macaque Genome Mmul_1"
    IRGSP_1_0 = "Ensembl Rice Genome IRGSP-1.0"
    CHIMP2_1_4 = "Ensembl Chimpanzee Genome CHIMP2.1.4"
    Rnor_5_0 = "Ensembl Rat Genome Rnor_5.0"
    Rnor_6_0 = "Ensembl Rat Genome Rnor_6.0"
    R64_1_1 = "Ensembl Yeast Genome R64-1-1"
    EF2 = "Ensembl Schizosaccharomyces pombe Genome EF2"
    Sbi1 = "Ensembl Sorghum Genome Sbi1"
    Sscrofa10_2 = "Ensembl Pig Genome Sscrofa10.2"
    AGPv3 = "Ensembl Maize Genome AGPv3"
    hg38 = "UCSC Human Genome hg38"
    hg19 = "UCSC Human Genome hg19"
    mm10 = "UCSC Mouse Genome mm10"
    bosTau8 = "UCSC Cow Genome bosTau8"
    ce10 = "UCSC Worm (C. elegans) Genome ce10"
    canFam3 = "UCSC Dog Genome canFam3"
    danRer10 = "UCSC Zebrafish Genome danRer10"
    dm6 = "UCSC Drosophila Genome dm6"
    equCab2 = "UCSC Horse Genome equCab2"
    galGal4 = "UCSC Chicken Genome galGal4"
    panTro4 = "UCSC Chimpanzee Genome panTro4"
    rn6 = "UCSC Rat Genome rn6"
    sacCer3 = "UCSC Yeast Genome sacCer3"
    susScr3 = "UCSC Pig Genome susScr3"


class StepOptions(Enum):
    mapping = "mapping"
    markduplicates = "markduplicates"
    prepare_recalibration = "prepare_recalibration"
    recalibrate = "recalibrate"
    variant_calling = "variant_calling"
    annotate = "annotate"


class ReferenceType(Enum):
    homo_sapiens = "Homo sapiens (RefSeq GRCh38.p14)"
    mus_musculus = "Mus musculus (RefSeq GRCm39)"
    rattus_norvegicus = "Rattus norvegicus (RefSeq GRCr8)"


flow = [
    Section(
        "Input",
        Params(
            "input",
            "step",
        ),
    ),
    Section(
        "Reference Genome",
        Fork(
            "genome_source",
            "",
            # latch_genome_source=ForkBranch(
            #     "Reference Genome",
            #     Params(
            #         "latch_genome",
            #     ),
            # ),
            igenomes_reference=ForkBranch(
                "iGenomes Reference",
                Params(
                    "genome",
                ),
                Spoiler(
                    "iGenomes Reference Options",
                    Params(
                        "igenomes_base",
                        "igenomes_ignore",
                    ),
                ),
            ),
            custom=ForkBranch(
                "Custom Reference Genome",
                Params(
                    "fasta",
                    "fasta_fai",
                ),
                Spoiler(
                    "Reference Genome Options",
                    Params(
                        "dbsnp_vqsr",
                        "known_indels_vqsr",
                        "known_snps",
                        "known_snps_tbi",
                        "known_snps_vqsr",
                        "ngscheckmate_bed",
                        "snpeff_db",
                        "snpeff_genome",
                        "vep_genome",
                        "vep_species",
                        "vep_cache_version",
                        "save_reference",
                        "build_only_index",
                        "download_cache",
                        "vep_cache",
                        "snpeff_cache",
                    ),
                ),
            ),
        ),
    ),
    Section(
        "Run Options",
        Params(
            "split_fastq",
            "wes",
            "intervals",
            "nucleotides_per_second",
            "no_intervals",
            "tools",
            "skip_tools",
        ),
    ),
    Section(
        "Output Directory",
        Params("run_name"),
        Text("Parent directory for outputs"),
        Params("outdir"),
    ),
    Spoiler(
        "Optional Arguments",
        Text("Additional optional arguments"),
        Spoiler(
            "FASTQ Preprocessing",
            Params(
                "trim_fastq",
                "umi_read_structure",
            ),
        ),
        Spoiler(
            "Preprocessing",
            Params(
                "aligner",
                "save_mapped",
                "save_output_as_bam",
                "use_gatk_spark",
            ),
        ),
        Spoiler(
            "Variant Calling",
            Params(
                "concatenate_vcfs",
                "only_paired_variant_calling",
                "joint_germline",
                "joint_mutect2",
            ),
        ),
        Spoiler(
            "Annotation",
            Params(
                "vep_custom_args",
                "vep_version",
                "bcftools_annotations",
                "bcftools_annotations_tbi",
                "bcftools_header_lines",
            ),
        ),
        Spoiler(
            "Reporting Options",
            Params(
                "email",
                "multiqc_title",
                "multiqc_methods_description",
            ),
        ),
    ),
]


generated_parameters = {
    "run_name": NextflowParameter(
        type=str,
        display_name="Run Name",
        description="Name of run",
        batch_table_column=True,
        rules=[
            LatchRule(
                regex=r"^[a-zA-Z0-9_-]+$",
                message="Run name must contain only letters, digits, underscores, and dashes. No spaces are allowed.",
            )
        ],
    ),
    "input": NextflowParameter(
        type=List[Sample],
        display_name="Input Samples",
        default=None,
        samplesheet=True,
        samplesheet_type="csv",
    ),
    "step": NextflowParameter(
        type=StepOptions,
        display_name="Starting Step",
        default="mapping",
        description="Starting step",
    ),
    "genome_source": NextflowParameter(
        type=str,
        display_name="Reference Genome",
        description="Choose Reference Genome",
    ),
    # "latch_genome": NextflowParameter(
    #     type=ReferenceType,
    #     display_name="Latch Verfied Reference Genome",
    #     description="Name of Latch Verfied Reference Genome.",
    #     default=ReferenceType.homo_sapiens,
    # ),
    "outdir": NextflowParameter(
        type=LatchOutputDir,
        display_name="Output Directory",
        default=None,
        description="The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.",
    ),
    "split_fastq": NextflowParameter(
        type=Optional[int],
        display_name="Split FastQ",
        default=50000000,
        description="Specify how many reads each split of a FastQ file contains. Set 0 to turn off splitting at all.",
    ),
    "wes": NextflowParameter(
        type=bool,
        display_name="Whole Exome Sequencing",
        default=None,
        description="Enable when exome or panel data is provided.",
    ),
    "intervals": NextflowParameter(
        type=Optional[LatchFile],
        display_name="Intervals File",
        default=None,
        description="Path to target bed file in case of whole exome or targeted sequencing or intervals file.",
    ),
    "nucleotides_per_second": NextflowParameter(
        type=Optional[int],
        display_name="Nucleotides per Second",
        default=200000,
        description="Estimate interval size.",
    ),
    "no_intervals": NextflowParameter(
        type=bool,
        display_name="Disable Intervals",
        default=None,
        description="Disable usage of intervals.",
    ),
    "tools": NextflowParameter(
        type=Optional[str],
        display_name="Tools to Use",
        default=None,
        description="Tools to use for duplicate marking, variant calling and/or for annotation.",
    ),
    "skip_tools": NextflowParameter(
        type=Optional[str],
        display_name="Tools to Skip",
        default=None,
        description="Disable specified tools.",
    ),
    "trim_fastq": NextflowParameter(
        type=bool,
        display_name="Trim FastQ",
        default=None,
        description="Run FastP for read trimming",
    ),
    "umi_read_structure": NextflowParameter(
        type=Optional[str],
        display_name="UMI Read Structure",
        default=None,
        description="Specify UMI read structure",
    ),
    "aligner": NextflowParameter(
        type=Optional[str],
        display_name="Aligner",
        default="bwa-mem",
        description="Specify aligner to be used to map reads to reference genome.",
    ),
    "save_mapped": NextflowParameter(
        type=bool,
        display_name="Save Mapped Files",
        default=None,
        description="Save mapped files.",
    ),
    "save_output_as_bam": NextflowParameter(
        type=bool,
        display_name="Save Output as BAM",
        default=None,
        description="Saves output from mapping (if `--save_mapped`), Markduplicates & Baserecalibration as BAM file instead of CRAM",
    ),
    "use_gatk_spark": NextflowParameter(
        type=Optional[str],
        display_name="Use GATK Spark",
        default=None,
        description="Enable usage of GATK Spark implementation for duplicate marking and/or base quality score recalibration",
    ),
    "concatenate_vcfs": NextflowParameter(
        type=bool,
        display_name="Concatenate VCFs",
        default=None,
        description="Option for concatenating germline vcf-files.",
    ),
    "only_paired_variant_calling": NextflowParameter(
        type=bool,
        display_name="Only Paired Variant Calling",
        default=None,
        description="If true, skips germline variant calling for matched normal to tumor sample. Normal samples without matched tumor will still be processed through germline variant calling tools.",
    ),
    "joint_germline": NextflowParameter(
        type=bool,
        display_name="Joint Germline Calling",
        default=None,
        description="Turn on the joint germline variant calling for GATK haplotypecaller",
    ),
    "joint_mutect2": NextflowParameter(
        type=bool,
        display_name="Joint Mutect2",
        default=None,
        description="Runs Mutect2 in joint (multi-sample) mode for better concordance among variant calls of tumor samples from the same patient. Mutect2 outputs will be stored in a subfolder named with patient ID under `variant_calling/mutect2/` folder. Only a single normal sample per patient is allowed. Tumor-only mode is also supported.",
    ),
    "vep_custom_args": NextflowParameter(
        type=Optional[str],
        display_name="VEP Custom Arguments",
        default="--everything --filter_common --per_gene --total_length --offline --format vcf",
        description="Add an extra custom argument to VEP.",
    ),
    "vep_version": NextflowParameter(
        type=Optional[str],
        display_name="VEP Version",
        default="111.0-0",
        description="Should reflect the VEP version used in the container.",
    ),
    "bcftools_annotations": NextflowParameter(
        type=Optional[LatchFile],
        display_name="BCFtools Annotations",
        default=None,
        description="A vcf file containing custom annotations to be used with bcftools annotate. Needs to be bgzipped.",
    ),
    "bcftools_annotations_tbi": NextflowParameter(
        type=Optional[LatchFile],
        display_name="BCFtools Annotations Index",
        default=None,
        description="Index file for `bcftools_annotations`",
    ),
    "bcftools_header_lines": NextflowParameter(
        type=Optional[LatchFile],
        display_name="BCFtools Header Lines",
        default=None,
        description="Text file with the header lines of `bcftools_annotations`",
    ),
    "genome": NextflowParameter(
        type=Optional[str],
        display_name="Reference Genome",
        description="Name of iGenomes reference.",
    ),
    "dbsnp_vqsr": NextflowParameter(
        type=Optional[str],
        display_name="dbSNP VQSR",
        default=None,
        description="label string for VariantRecalibration (haplotypecaller joint variant calling)",
    ),
    "fasta": NextflowParameter(
        type=Optional[LatchFile],
        display_name="FASTA Genome File",
        default=None,
        description="Path to FASTA genome file.",
    ),
    "fasta_fai": NextflowParameter(
        type=Optional[str],
        display_name="FASTA Index",
        default=None,
        description="Path to FASTA reference index.",
    ),
    "known_indels_vqsr": NextflowParameter(
        type=Optional[str],
        display_name="Known Indels VQSR",
        default=None,
        description="1st label string for VariantRecalibration (haplotypecaller joint variant calling)",
    ),
    "known_snps": NextflowParameter(
        type=Optional[LatchFile],
        display_name="Known SNPs",
        default=None,
        description="Path to known snps file.",
    ),
    "known_snps_tbi": NextflowParameter(
        type=Optional[LatchFile],
        display_name="Known SNPs Index",
        default=None,
        description="Path to known snps file index.",
    ),
    "known_snps_vqsr": NextflowParameter(
        type=Optional[str],
        display_name="Known SNPs VQSR",
        default=None,
        description="label string for VariantRecalibration (haplotypecaller joint variant calling)",
    ),
    "ngscheckmate_bed": NextflowParameter(
        type=Optional[LatchFile],
        display_name="NGSCheckMate BED",
        default=None,
        description="Path to SNP bed file for sample checking with NGSCheckMate",
    ),
    "snpeff_db": NextflowParameter(
        type=Optional[str],
        display_name="snpEff DB Version",
        default=None,
        description="snpEff DB version.",
    ),
    "snpeff_genome": NextflowParameter(
        type=Optional[str],
        display_name="snpEff Genome",
        default=None,
        description="snpEff genome.",
    ),
    "vep_genome": NextflowParameter(
        type=Optional[str],
        display_name="VEP Genome",
        default=None,
        description="VEP genome.",
    ),
    "vep_species": NextflowParameter(
        type=Optional[str],
        display_name="VEP Species",
        default=None,
        description="VEP species.",
    ),
    "vep_cache_version": NextflowParameter(
        type=Optional[str],
        display_name="VEP Cache Version",
        default=None,
        description="VEP cache version.",
    ),
    "save_reference": NextflowParameter(
        type=bool,
        display_name="Save References",
        default=None,
        description="Save built references.",
    ),
    "build_only_index": NextflowParameter(
        type=bool,
        display_name="Build Only Index",
        default=None,
        description="Only built references.",
    ),
    "download_cache": NextflowParameter(
        type=bool,
        display_name="Download Cache",
        default=None,
        description="Download annotation cache.",
    ),
    "igenomes_base": NextflowParameter(
        type=Optional[LatchDir],
        display_name="iGenomes Base",
        default=None,
        description="Directory / URL base for iGenomes references.",
    ),
    "igenomes_ignore": NextflowParameter(
        type=bool,
        display_name="Ignore iGenomes",
        default=None,
        description="Do not load the iGenomes reference config.",
    ),
    "vep_cache": NextflowParameter(
        type=Optional[LatchDir],
        display_name="VEP Cache",
        default=None,
        description="Path to VEP cache.",
    ),
    "snpeff_cache": NextflowParameter(
        type=Optional[LatchDir],
        display_name="snpEff Cache",
        default=None,
        description="Path to snpEff cache.",
    ),
    "email": NextflowParameter(
        type=Optional[str],
        display_name="Email",
        default=None,
        description="Email address for completion summary.",
    ),
    "multiqc_title": NextflowParameter(
        type=Optional[str],
        display_name="MultiQC Title",
        default=None,
        description="MultiQC report title. Printed as page header, used for filename if not otherwise specified.",
    ),
    "multiqc_methods_description": NextflowParameter(
        type=Optional[str],
        display_name="MultiQC Methods Description",
        default=None,
        description="Custom MultiQC yaml file containing HTML including a methods description.",
    ),
}
