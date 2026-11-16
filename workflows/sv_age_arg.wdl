version 1.0

workflow SVAgeARG {
    meta {
        description: "Dates one structural variant from a local ARG: tsinfer on the phased SNVs flanking the SV (its own span excluded), singletons added for dating, tsdate, then parsimony mapping of the SV onto the trees with genotype/phasing/recurrence QC."
    }

    parameter_meta {
        SvIds: "VCF entries that together make up the SV (an MEI split across several records is one SV): each is the VCF ID, CHROM:POS:REF:ALT or CHROM:POS_REF_ALT. A haplotype carries the SV if it carries any of them. Looked up in the phased VCF first, then the full VCF."
        Chrom: "Chromosome, named as in the VCFs."
        SvStart: "1-based start of the SV. SNVs between SvStart - SvExcludePadBp and SvEnd + SvExcludePadBp are never used, so a deletion's hemizygous span cannot distort the trees."
        SvEnd: "1-based end of the SV (equal to SvStart for an insertion)."
        PhasedVcfList: "Two-column table (tab or space separated, header optional): chromosome, path of the phased MAC >= 2 VCF for that chromosome. Chromosomes may be written chr7 or 7. Only the row for Chrom is used; that VCF is read remotely by region (never localized), so each must be bgzipped and indexed."
        FullVcfList: "Same layout, for the full callset that includes singletons (phasing not required)."
        PhasedVcfIndexList: "Optional table of the same layout giving index paths, for indexes not at <vcf>.tbi / <vcf>.csi."
        FullVcfIndexList: "Optional index table for the full VCFs."
        Prefix: "Output prefix; defaults to the first SvId with unsafe characters replaced."
        FlankBp: "SNVs this far either side of the SV are used to infer the ARG."
        SvExcludePadBp: "Extra exclusion around the SV span, covering imprecise breakpoints."
        SvSearchPadBp: "How far beyond SvStart/SvEnd to look for the SvIds records."
        AncestralFastaList: "Two-column table like PhasedVcfList: chromosome, path (gs://, https:// or local) of that chromosome's ancestral-allele FASTA, e.g. Ensembl's homo_sapiens_ancestor_7.fa. Must be indexed (<fasta>.fai beside it, plus .gzi if bgzipped; plain gzip cannot be read by region). Only the flank region of Chrom's file is read, remotely, never the whole file. Ensembl ANCESTOR_for_chromosome headers are recognised. Without it (or AncestralInfoField) REF is treated as ancestral and the result is flagged."
        AncestralInfoField: "INFO field holding the ancestral allele (e.g. AA); used only if AncestralFastaList is not given."
        AncestralHighConfidenceOnly: "Drop SNVs whose ancestral call is lower-case (low confidence) in the FASTA."
        SvType: "INS or DEL. INS: REF is the ancestral allele and the insertion (ALT) is dated. DEL: ALT is the ancestral allele, i.e. the reference sequence (e.g. a reference SVA missing in some people) is the derived insertion and is what gets dated."
        SvAncestralState: "Optional override of SvType's polarity: ref, alt, or auto (whichever parsimony prefers)."
        GeneticMap: "Optional genetic map: chrom/pos/rate/cM (HapMap/SHAPEIT style, header optional), PLINK .map, or pos/rate/cM; may hold all chromosomes. Without it RecombinationRate is used."
        MutationRate: "Per-bp per-generation rate used by tsdate."
        GenerationTime: "Years per generation for the *_years columns."
        SingletonsPhased: "Set true only if the full VCF's singleton calls are phased; otherwise tsdate treats singleton phase as unknown."
        ScanBp: "Also map the SV on every tree within this distance, as a robustness check."
        BackgroundBp: "SNVs of similar derived count between ScanBp and this distance are held out of inference and mapped afterwards like the SV, giving the background parsimony cost and clade-split rate the SV is compared against."
        MaxErrorCladeHaps: "An extra mutation explaining at most this many haplotypes is reported as a genotyping/phasing error."
        RecurrentMinHaps: "An extra forward mutation carried by at least this many haplotypes is reported as a possible recurrent event."
        RequesterPaysProject: "Billing project, if the VCF buckets are requester-pays."
        ImageTag: "Tag of ayenkin1871/mei-lr-association-arg (latest or a 7-char commit SHA)."
    }

    input {
        Array[String] SvIds
        String Chrom
        Int SvStart
        Int SvEnd
        File PhasedVcfList
        File FullVcfList
        File? PhasedVcfIndexList
        File? FullVcfIndexList
        String? Prefix

        Int FlankBp = 500000
        Int SvExcludePadBp = 500
        Int SvSearchPadBp = 1000

        File? AncestralFastaList
        String? AncestralInfoField
        Boolean AncestralHighConfidenceOnly = false
        String SvType = "INS"
        String? SvAncestralState

        File? GeneticMap
        Float RecombinationRate = 0.00000001
        Float MismatchRatio = 1.0
        Float MutationRate = 0.0000000129
        Float GenerationTime = 28.0
        Boolean SingletonsPhased = false

        Int ScanBp = 20000
        Int BackgroundBp = 100000
        Int MaxErrorCladeHaps = 2
        Int RecurrentMinHaps = 5

        String? RequesterPaysProject
        String ImageTag = "latest"
        Int Cpu = 8
        # 20k haplotypes x 14k SNVs over 1 Mb peaked at 1.1 GB and ran 16 min on 8 threads
        Int MemoryGB = 8
        Int DiskGB = 50
        Int Preemptible = 1
    }

    String prefix = select_first([Prefix, sub(SvIds[0], "[^A-Za-z0-9._-]", "_")])
    String docker = "ayenkin1871/mei-lr-association-arg:" + ImageTag

    call ExtractRegion {
        input:
            PhasedVcfList = PhasedVcfList,
            FullVcfList = FullVcfList,
            PhasedVcfIndexList = PhasedVcfIndexList,
            FullVcfIndexList = FullVcfIndexList,
            AncestralFastaList = AncestralFastaList,
            Chrom = Chrom,
            SvStart = SvStart,
            SvEnd = SvEnd,
            FlankBp = FlankBp,
            SvSearchPadBp = SvSearchPadBp,
            RequesterPaysProject = RequesterPaysProject,
            docker = docker,
            DiskGB = DiskGB,
            Preemptible = Preemptible
    }

    call InferAndDate {
        input:
            RegionPhasedVcf = ExtractRegion.RegionPhasedVcf,
            RegionSingletonVcf = ExtractRegion.RegionSingletonVcf,
            SvPhasedVcf = ExtractRegion.SvPhasedVcf,
            SvFullVcf = ExtractRegion.SvFullVcf,
            SvIds = SvIds,
            Chrom = Chrom,
            SvStart = SvStart,
            SvEnd = SvEnd,
            FlankBp = FlankBp,
            SvExcludePadBp = SvExcludePadBp,
            AncestralFasta = ExtractRegion.AncestralFasta,
            AncestralOffset = ExtractRegion.AncestralOffset,
            AncestralInfoField = AncestralInfoField,
            AncestralHighConfidenceOnly = AncestralHighConfidenceOnly,
            SvType = SvType,
            SvAncestralState = SvAncestralState,
            GeneticMap = GeneticMap,
            RecombinationRate = RecombinationRate,
            MismatchRatio = MismatchRatio,
            MutationRate = MutationRate,
            GenerationTime = GenerationTime,
            SingletonsPhased = SingletonsPhased,
            ScanBp = ScanBp,
            BackgroundBp = BackgroundBp,
            MaxErrorCladeHaps = MaxErrorCladeHaps,
            RecurrentMinHaps = RecurrentMinHaps,
            Prefix = prefix,
            docker = docker,
            Cpu = Cpu,
            MemoryGB = MemoryGB,
            DiskGB = DiskGB,
            Preemptible = Preemptible
    }

    output {
        File Arg = InferAndDate.Arg
        File SvAge = InferAndDate.SvAge
        File SvMutations = InferAndDate.SvMutations
        File SuspiciousSamples = InferAndDate.SuspiciousSamples
        File SvGenotypes = InferAndDate.SvGenotypes
        File TreeScan = InferAndDate.TreeScan
        File SvEntries = InferAndDate.SvEntries
        File RegionQc = InferAndDate.RegionQc
        File Log = InferAndDate.Log
        File ExtractCounts = ExtractRegion.Counts
    }
}

task ExtractRegion {
    input {
        File PhasedVcfList
        File FullVcfList
        File? PhasedVcfIndexList
        File? FullVcfIndexList
        File? AncestralFastaList
        String Chrom
        Int SvStart
        Int SvEnd
        Int FlankBp
        Int SvSearchPadBp
        String? RequesterPaysProject
        String docker
        Int DiskGB
        Int Preemptible
    }

    command <<<
        set -euo pipefail

        chrom="~{Chrom}"
        chrom_n=${chrom#chr}
        lookup() {
            local list=$1 required=$2 v
            [[ -z "$list" ]] && return 0
            v=$(tr -d '\r' < "$list" | awk -v a="$chrom" -v b="$chrom_n" -v c="chr$chrom_n" \
                '$1==a{print $2; f=1; exit} $1==b||$1==c{x=$2} END{if(!f && x!="") print x}')
            if [[ -z "$v" && "$required" == 1 ]]; then
                echo "No row for $chrom in $list; chromosomes listed: $(tr -d '\r' < "$list" | awk 'NF{print $1}' | tr '\n' ' ')" >&2
                exit 1
            fi
            printf '%s' "$v"
        }
        phased=$(lookup "~{PhasedVcfList}" 1)
        full=$(lookup "~{FullVcfList}" 1)
        phased_idx=$(lookup "~{default="" PhasedVcfIndexList}" 0)
        full_idx=$(lookup "~{default="" FullVcfIndexList}" 0)
        [[ -n "$phased_idx" ]] && phased="${phased}##idx##${phased_idx}"
        [[ -n "$full_idx" ]] && full="${full}##idx##${full_idx}"
        anc=$(lookup "~{default="" AncestralFastaList}" 1)
        echo "phased VCF: $phased"
        echo "full VCF:   $full"
        echo "ancestral:  ${anc:-none}"

        if [[ "$phased$full$anc" == *gs://* ]]; then
            token=$(curl -sf -H "Metadata-Flavor: Google" \
                "http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token" \
                | python3 -c 'import json,sys; print(json.load(sys.stdin)["access_token"])') \
                || { echo "Could not get a GCS access token from the metadata server" >&2; exit 1; }
            export GCS_OAUTH_TOKEN="$token"
        fi
        ~{if defined(RequesterPaysProject) then "export GCS_REQUESTER_PAYS_PROJECT=" + RequesterPaysProject else ""}


        lo=$(( ~{SvStart} - ~{FlankBp} )); (( lo < 1 )) && lo=1
        hi=$(( ~{SvEnd} + ~{FlankBp} ))
        slo=$(( ~{SvStart} - ~{SvSearchPadBp} )); (( slo < 1 )) && slo=1
        shi=$(( ~{SvEnd} + ~{SvSearchPadBp} ))

        # only the flank region of the ancestral FASTA is read, remotely like the VCFs; the
        # .fai is fetched first because Ensembl contig names (ANCESTOR_for_chromosome:GRCh38:7:...)
        # contain colons and must be named exactly in {name}:start-end form
        : > ancestral.region.fa
        echo 0 > ancestral_offset.txt
        if [[ -n "$anc" ]]; then
            case "$anc" in
                gs://*)
                    hdr=(-H "Authorization: Bearer $GCS_OAUTH_TOKEN")
                    [[ -n "${GCS_REQUESTER_PAYS_PROJECT:-}" ]] && \
                        hdr+=(-H "x-goog-user-project: $GCS_REQUESTER_PAYS_PROJECT")
                    curl -sSf "${hdr[@]}" -o ancestral.fai \
                        "https://storage.googleapis.com/${anc#gs://}.fai" ;;
                http://*|https://*) curl -sSfL -o ancestral.fai "$anc.fai" ;;
                *) cp "$anc.fai" ancestral.fai ;;
            esac || { echo "Could not read $anc.fai: the ancestral FASTA must be indexed (samtools faidx)" >&2; exit 1; }
            contig=$(awk -F'\t' -v a="$chrom" -v b="$chrom_n" '
                { split($1, p, ":") }
                $1 == a || $1 == b || $1 == "chr" b || (p[1] ~ /^ANCESTOR/ && (p[3] == b || "chr" p[3] == a)) { print $1; exit }
                ' ancestral.fai)
            if [[ -z "$contig" && $(grep -c . ancestral.fai) -eq 1 ]]; then
                contig=$(cut -f1 ancestral.fai)
            fi
            if [[ -z "$contig" ]]; then
                echo "No contig for $chrom in $anc.fai; contigs: $(cut -f1 ancestral.fai | tr '\n' ' ')" >&2
                exit 1
            fi
            samtools faidx "$anc" "{$contig}:$lo-$hi" > ancestral.region.fa
            echo "$lo" > ancestral_offset.txt
            echo "ancestral: $contig:$lo-$hi ($(grep -v '^>' ancestral.region.fa | tr -d '\n' | wc -c) bp)"
        fi

        bcftools query -l "$phased" > phased_samples.txt
        bcftools view --no-version -r "$chrom:$lo-$hi" -v snps -m2 -M2 "$phased" \
            -Oz -o region.phased.vcf.gz
        # -C 1:minor is evaluated after -S, i.e. among the phased VCF's samples only
        bcftools view --no-version -r "$chrom:$lo-$hi" -v snps -m2 -M2 \
            -S phased_samples.txt --force-samples -C 1:minor "$full" \
            -Oz -o region.singletons.vcf.gz
        bcftools view --no-version -r "$chrom:$slo-$shi" "$phased" -Oz -o sv.phased.vcf.gz
        bcftools view --no-version -r "$chrom:$slo-$shi" "$full" -Oz -o sv.full.vcf.gz

        printf 'file\trecords\n' > extract_counts.tsv
        for f in region.phased region.singletons sv.phased sv.full; do
            printf '%s\t%s\n' "$f" "$(bcftools view -H "$f.vcf.gz" | wc -l)" >> extract_counts.tsv
        done
        cat extract_counts.tsv
        if [[ $(awk '$1=="region.phased"{print $2}' extract_counts.tsv) -eq 0 ]]; then
            echo "No biallelic SNVs in $chrom:$lo-$hi of the phased VCF (check Chrom naming and the template)" >&2
            exit 1
        fi
    >>>

    runtime {
        docker: docker
        memory: "4G"
        cpu: 2
        disks: "local-disk " + DiskGB + " SSD"
        preemptible: Preemptible
        maxRetries: 1
    }

    output {
        File RegionPhasedVcf = "region.phased.vcf.gz"
        File RegionSingletonVcf = "region.singletons.vcf.gz"
        File SvPhasedVcf = "sv.phased.vcf.gz"
        File SvFullVcf = "sv.full.vcf.gz"
        File Counts = "extract_counts.tsv"
        File AncestralFasta = "ancestral.region.fa"
        Int AncestralOffset = read_int("ancestral_offset.txt")
    }
}

task InferAndDate {
    input {
        File RegionPhasedVcf
        File RegionSingletonVcf
        File SvPhasedVcf
        File SvFullVcf
        Array[String] SvIds
        String Chrom
        Int SvStart
        Int SvEnd
        Int FlankBp
        Int SvExcludePadBp
        File AncestralFasta
        Int AncestralOffset
        String? AncestralInfoField
        Boolean AncestralHighConfidenceOnly
        String SvType
        String? SvAncestralState
        File? GeneticMap
        Float RecombinationRate
        Float MismatchRatio
        Float MutationRate
        Float GenerationTime
        Boolean SingletonsPhased
        Int ScanBp
        Int BackgroundBp
        Int MaxErrorCladeHaps
        Int RecurrentMinHaps
        String Prefix
        String docker
        Int Cpu
        Int MemoryGB
        Int DiskGB
        Int Preemptible
    }

    command <<<
        set -euo pipefail
        anc_args=()
        if [[ -s "~{AncestralFasta}" ]]; then
            anc_args=(--ancestral-fasta "~{AncestralFasta}" --ancestral-offset ~{AncestralOffset})
        fi
        python3 /scripts/sv_age_arg.py \
            --phased-vcf ~{RegionPhasedVcf} \
            --singleton-vcf ~{RegionSingletonVcf} \
            --sv-phased-vcf ~{SvPhasedVcf} \
            --sv-full-vcf ~{SvFullVcf} \
            --sv-ids '~{sep="," SvIds}' \
            --chrom ~{Chrom} \
            --sv-start ~{SvStart} \
            --sv-end ~{SvEnd} \
            --flank ~{FlankBp} \
            --exclude-pad ~{SvExcludePadBp} \
            "${anc_args[@]}" \
            ~{"--ancestral-info-field " + AncestralInfoField} \
            ~{if AncestralHighConfidenceOnly then "--ancestral-high-confidence-only" else ""} \
            --sv-type ~{SvType} \
            ~{"--sv-ancestral " + SvAncestralState} \
            ~{"--genetic-map " + GeneticMap} \
            --recombination-rate ~{RecombinationRate} \
            --mismatch-ratio ~{MismatchRatio} \
            --mutation-rate ~{MutationRate} \
            --generation-time ~{GenerationTime} \
            ~{if SingletonsPhased then "--singletons-phased" else ""} \
            --scan-bp ~{ScanBp} \
            --background-bp ~{BackgroundBp} \
            --max-error-clade-haps ~{MaxErrorCladeHaps} \
            --recurrent-min-haps ~{RecurrentMinHaps} \
            --threads ~{Cpu} \
            --prefix ~{Prefix}
    >>>

    runtime {
        docker: docker
        memory: MemoryGB + "G"
        cpu: Cpu
        disks: "local-disk " + DiskGB + " SSD"
        preemptible: Preemptible
        maxRetries: 1
    }

    output {
        File Arg = Prefix + ".arg.trees.tsz"
        File SvAge = Prefix + ".sv_age.tsv"
        File SvMutations = Prefix + ".sv_mutations.tsv"
        File SuspiciousSamples = Prefix + ".suspicious_samples.tsv"
        File SvGenotypes = Prefix + ".sv_genotypes.tsv.gz"
        File TreeScan = Prefix + ".tree_scan.tsv"
        File SvEntries = Prefix + ".sv_entries.tsv"
        File RegionQc = Prefix + ".region_qc.tsv"
        File Log = Prefix + ".log"
    }
}
