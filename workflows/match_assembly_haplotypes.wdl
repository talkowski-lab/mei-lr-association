version 1.0

## Matches an assembly's two haplotype BAMs (aligned to the reference, e.g. by
## AlignDiploidAssemblyToRef) to the phasing in a phased VCF for the same
## sample over a region of interest. For every phased, heterozygous, biallelic
## SNV the sample has in the region, reads the base each assembly haplotype
## carries at that position (no assembly variant calling needed) and reports
## percent concordance of each assembly haplotype against each VCF-phased
## haplotype, per phase set.

workflow MatchAssemblyHaplotypes {
    input {
        File Hap1Bam
        File Hap1BamIndex
        File Hap2Bam
        File Hap2BamIndex
        File Vcf
        File VcfIndex
        String SampleId
        String Region
        String Prefix
        String ImageTag = "latest"
        Int MinMapQ = 5
    }

    call MatchHaplotypes {
        input:
            Hap1Bam = Hap1Bam,
            Hap1BamIndex = Hap1BamIndex,
            Hap2Bam = Hap2Bam,
            Hap2BamIndex = Hap2BamIndex,
            Vcf = Vcf,
            VcfIndex = VcfIndex,
            SampleId = SampleId,
            Region = Region,
            Prefix = Prefix,
            ImageTag = ImageTag,
            MinMapQ = MinMapQ
    }

    output {
        File PerSiteTable = MatchHaplotypes.PerSiteTable
        File SummaryTable = MatchHaplotypes.SummaryTable
    }
}

task MatchHaplotypes {
    input {
        File Hap1Bam
        File Hap1BamIndex
        File Hap2Bam
        File Hap2BamIndex
        File Vcf
        File VcfIndex
        String SampleId
        String Region
        String Prefix
        String ImageTag = "latest"
        Int MinMapQ = 5
    }

    Int disk_size = ceil(size(Hap1Bam, "GB") + size(Hap2Bam, "GB") + size(Vcf, "GB") * 1.5) + 20

    command <<<
        set -euo pipefail

        # Co-locate indexes next to their files so pysam/tabix find them.
        ln -s ~{Hap1Bam} hap1.bam
        ln -s ~{Hap1BamIndex} hap1.bam.bai
        ln -s ~{Hap2Bam} hap2.bam
        ln -s ~{Hap2BamIndex} hap2.bam.bai
        ln -s ~{Vcf} input.vcf.gz
        ln -s ~{VcfIndex} input.vcf.gz.tbi

        python3 /scripts/match_assembly_haplotypes.py \
            --hap1-bam hap1.bam \
            --hap2-bam hap2.bam \
            --vcf input.vcf.gz \
            --sample-id ~{SampleId} \
            --region ~{Region} \
            --prefix ~{Prefix} \
            --min-mapq ~{MinMapQ}
    >>>

    runtime {
        docker: "ayenkin1871/mei-lr-association-python_general:" + ImageTag
        memory: "8G"
        cpu: 2
        disks: "local-disk " + disk_size + " SSD"
        preemptible: 3
        maxRetries: 2
    }

    output {
        File PerSiteTable = "~{Prefix}.persite.tsv"
        File SummaryTable = "~{Prefix}.summary.tsv"
    }
}
