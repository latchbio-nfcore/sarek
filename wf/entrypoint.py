import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import List, Optional

import requests
from latch.executions import rename_current_execution, report_nextflow_used_storage
from latch.ldata.path import LPath
from latch.resources.tasks import custom_task, nextflow_runtime_task
from latch.types import metadata
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile
from latch_cli.nextflow.utils import _get_execution_name
from latch_cli.nextflow.workflow import get_flag
from latch_cli.services.register.utils import import_module_by_path
from latch_cli.utils import urljoins

meta = Path("latch_metadata") / "__init__.py"
import_module_by_path(meta)

sys.stdout.reconfigure(line_buffering=True)


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


class Aligner(Enum):
    bwa_mem = "bwa-mem"
    bwa_mem2 = "bwa-mem2"
    dragmap = "dragmap"
    sentieon_bwamen = "sentieon-bwamem"


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


@custom_task(cpu=0.25, memory=0.5, storage_gib=1)
def initialize(run_name: str) -> str:
    rename_current_execution(str(run_name))

    token = os.environ.get("FLYTE_INTERNAL_EXECUTION_ID")
    if token is None:
        raise RuntimeError("failed to get execution token")

    headers = {"Authorization": f"Latch-Execution-Token {token}"}

    print("Provisioning shared storage volume... ", end="")
    resp = requests.post(
        "http://nf-dispatcher-service.flyte.svc.cluster.local/provision-storage-ofs",
        # "http://nf-dispatcher-service.flyte.svc.cluster.local/provision-storage",
        headers=headers,
        json={
            "storage_expiration_hours": 120,
            "version": 2,
        },
    )
    resp.raise_for_status()
    print("Done.")

    return resp.json()["name"]


input_construct_samplesheet = metadata._nextflow_metadata.parameters[
    "input"
].samplesheet_constructor


@nextflow_runtime_task(cpu=4, memory=8, storage_gib=100)
def nextflow_runtime(
    pvc_name: str,
    run_name: str,
    input: List[Sample],
    genome_source: str,
    # latch_genome: ReferenceType,
    genome: Optional[GenomeReference],
    outdir: LatchOutputDir,
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
    fasta_fai: Optional[str],
    known_indels_vqsr: Optional[str],
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
    multiqc_methods_description: Optional[str],
    step: StepOptions,
    split_fastq: int,
    nucleotides_per_second: int,
    aligner: Aligner,
    vep_custom_args: Optional[str],
    vep_version: Optional[str],
) -> None:
    shared_dir = Path("/nf-workdir")

    input_samplesheet = input_construct_samplesheet(input)

    ignore_list = [
        "latch",
        ".latch",
        ".git",
        "nextflow",
        ".nextflow",
        "work",
        "results",
        "miniconda",
        "anaconda3",
        "mambaforge",
    ]

    shutil.copytree(
        Path("/root"),
        shared_dir,
        ignore=lambda src, names: ignore_list,
        ignore_dangling_symlinks=True,
        dirs_exist_ok=True,
    )

    cmd = [
        "/root/nextflow",
        "run",
        str(shared_dir / "main.nf"),
        "-work-dir",
        str(shared_dir),
        "-profile",
        "docker",
        "-c",
        "latch.config",
        "-resume",
        *get_flag("input", input_samplesheet),
        *get_flag("step", step),
        *get_flag("outdir", LatchOutputDir(f"{outdir.remote_path}/{run_name}")),
        *get_flag("split_fastq", split_fastq),
        *get_flag("wes", wes),
        *get_flag("intervals", intervals),
        *get_flag("nucleotides_per_second", nucleotides_per_second),
        *get_flag("no_intervals", no_intervals),
        *get_flag("tools", tools),
        *get_flag("skip_tools", skip_tools),
        *get_flag("trim_fastq", trim_fastq),
        *get_flag("umi_read_structure", umi_read_structure),
        *get_flag("aligner", aligner),
        *get_flag("save_mapped", save_mapped),
        *get_flag("save_output_as_bam", save_output_as_bam),
        *get_flag("use_gatk_spark", use_gatk_spark),
        *get_flag("concatenate_vcfs", concatenate_vcfs),
        *get_flag("only_paired_variant_calling", only_paired_variant_calling),
        *get_flag("joint_germline", joint_germline),
        *get_flag("joint_mutect2", joint_mutect2),
        *get_flag("vep_custom_args", vep_custom_args),
        *get_flag("vep_version", vep_version),
        *get_flag("bcftools_annotations", bcftools_annotations),
        *get_flag("bcftools_annotations_tbi", bcftools_annotations_tbi),
        *get_flag("bcftools_header_lines", bcftools_header_lines),
        *get_flag("genome", str(genome.name.replace("_", ".")) if genome else None),
        *get_flag("fasta", fasta),
        *get_flag("fasta_fai", fasta_fai),
        *get_flag("dbsnp_vqsr", dbsnp_vqsr),
        *get_flag("known_indels_vqsr", known_indels_vqsr),
        *get_flag("known_snps", known_snps),
        *get_flag("known_snps_tbi", known_snps_tbi),
        *get_flag("known_snps_vqsr", known_snps_vqsr),
        *get_flag("ngscheckmate_bed", ngscheckmate_bed),
        *get_flag("snpeff_db", snpeff_db),
        *get_flag("snpeff_genome", snpeff_genome),
        *get_flag("vep_genome", vep_genome),
        *get_flag("vep_species", vep_species),
        *get_flag("vep_cache_version", vep_cache_version),
        *get_flag("save_reference", save_reference),
        *get_flag("build_only_index", build_only_index),
        *get_flag("download_cache", download_cache),
        *get_flag("igenomes_base", igenomes_base),
        *get_flag("igenomes_ignore", igenomes_ignore),
        *get_flag("vep_cache", vep_cache),
        *get_flag("snpeff_cache", snpeff_cache),
        *get_flag("email", email),
        *get_flag("multiqc_title", multiqc_title),
        *get_flag("multiqc_methods_description", multiqc_methods_description),
    ]

    # if genome_source == "latch_genome_source":
    #     cmd += [
    #         "--fasta",
    #         f"s3://latch-public/nf-core/rnaseq/{latch_genome.name}/{latch_genome.name}.genomic.fna",
    #         "--bismark_index",
    #         f"s3://latch-public/nf-core/rnaseq/{latch_genome.name}/{latch_genome.name}.genomic.gtf",
    #     ]

    print("Launching Nextflow Runtime")
    print(" ".join(cmd))
    print(flush=True)

    failed = False
    try:
        env = {
            **os.environ,
            "NXF_ANSI_LOG": "false",
            "NXF_HOME": "/root/.nextflow",
            "NXF_OPTS": "-Xms1536M -Xmx6144M -XX:ActiveProcessorCount=4",
            "NXF_DISABLE_CHECK_LATEST": "true",
            "NXF_ENABLE_VIRTUAL_THREADS": "false",
        }
        subprocess.run(
            cmd,
            env=env,
            check=True,
            cwd=str(shared_dir),
        )
    except subprocess.CalledProcessError:
        failed = True
    finally:
        print()

        nextflow_log = shared_dir / ".nextflow.log"
        if nextflow_log.exists():
            name = _get_execution_name()
            if name is None:
                print("Skipping logs upload, failed to get execution name")
            else:
                remote = LPath(
                    urljoins(
                        "latch:///your_log_dir/nf_nf_core_sarek", name, "nextflow.log"
                    )
                )
                print(f"Uploading .nextflow.log to {remote.path}")
                remote.upload_from(nextflow_log)

        print("Computing size of workdir... ", end="")
        try:
            result = subprocess.run(
                ["du", "-sb", str(shared_dir)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=5 * 60,
            )

            size = int(result.stdout.split()[0])
            report_nextflow_used_storage(size)
            print(f"Done. Workdir size: {size / 1024 / 1024 / 1024: .2f} GiB")
        except subprocess.TimeoutExpired:
            print(
                "Failed to compute storage size: Operation timed out after 5 minutes."
            )
        except subprocess.CalledProcessError as e:
            print(f"Failed to compute storage size: {e.stderr}")
        except Exception as e:
            print(f"Failed to compute storage size: {e}")

    if failed:
        sys.exit(1)
