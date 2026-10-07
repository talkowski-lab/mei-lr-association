version 1.0

# Given a gzipped sparse long-form genotype table with columns including
# `indiv` and `var_id`, emit the subset of a variant BED whose name column is
# present in `var_id` for one requested individual.
workflow IndividualVariantBedFromLongGenotypes {
    input {
        File GenotypeTable
        File VariantsBed
        String Individual
        Int LocusCol = 4
        String Prefix = Individual
        String DockerImage = "debian:bookworm-slim"
        Int MemoryGB = 2
        Int CPU = 1
        Int? DiskGB
    }

    call FilterVariantBedByIndividual {
        input:
            GenotypeTable = GenotypeTable,
            VariantsBed = VariantsBed,
            Individual = Individual,
            LocusCol = LocusCol,
            Prefix = Prefix,
            DockerImage = DockerImage,
            MemoryGB = MemoryGB,
            CPU = CPU,
            DiskGB = DiskGB
    }

    output {
        File VariantBed = FilterVariantBedByIndividual.VariantBed
    }
}

task FilterVariantBedByIndividual {
    input {
        File GenotypeTable
        File VariantsBed
        String Individual
        Int LocusCol = 4
        String Prefix = Individual
        String DockerImage = "debian:bookworm-slim"
        Int MemoryGB = 2
        Int CPU = 1
        Int? DiskGB
    }

    Int auto_disk_size = ceil((size(GenotypeTable, "GB") + size(VariantsBed, "GB")) * 2) + 5

    command <<<
        set -euo pipefail

        if [[ ~{LocusCol} -lt 1 ]]; then
            echo "LocusCol must be >= 1, got ~{LocusCol}" >&2
            exit 1
        fi

        gzip -cd ~{GenotypeTable} > genotypes.tsv

        awk -F'\t' -v OFS='\t' -v indiv="~{Individual}" -v bed_name_col="~{LocusCol}" '
            NR == FNR {
                if (FNR == 1) {
                    indiv_col = 0
                    var_id_col = 0
                    for (i = 1; i <= NF; i++) {
                        if ($i == "indiv") {
                            indiv_col = i
                        } else if ($i == "var_id") {
                            var_id_col = i
                        }
                    }
                    if (indiv_col == 0 || var_id_col == 0) {
                        printf "GenotypeTable header must include indiv and var_id columns\n" > "/dev/stderr"
                        exit 1
                    }
                    next
                }
                if ($indiv_col == indiv && $var_id_col != "") {
                    variants[$var_id_col] = 1
                }
                next
            }
            /^#/ {
                print
                next
            }
            NF < bed_name_col {
                printf "VariantsBed row has %d columns, but LocusCol=%d: %s\n", NF, bed_name_col, $0 > "/dev/stderr"
                exit 1
            }
            $bed_name_col in variants {
                print
            }
        ' genotypes.tsv ~{VariantsBed} > "~{Prefix}.individual_variants.bed"

        n_var_ids=$(awk -F'\t' -v indiv="~{Individual}" '
            FNR == 1 {
                indiv_col = 0
                var_id_col = 0
                for (i = 1; i <= NF; i++) {
                    if ($i == "indiv") {
                        indiv_col = i
                    } else if ($i == "var_id") {
                        var_id_col = i
                    }
                }
                if (indiv_col == 0 || var_id_col == 0) {
                    exit 1
                }
                next
            }
            $indiv_col == indiv && $var_id_col != "" { seen[$var_id_col] = 1 }
            END { print length(seen) + 0 }
        ' genotypes.tsv)
        n_variants=$(awk 'BEGIN{n=0} !/^#/ {n++} END{print n}' "~{Prefix}.individual_variants.bed")

        echo "Found ${n_var_ids} unique non-reference var_id values for ~{Individual}" >&2
        echo "Wrote ${n_variants} matching BED rows to ~{Prefix}.individual_variants.bed" >&2
    >>>

    runtime {
        docker: DockerImage
        memory: MemoryGB + " GB"
        cpu: CPU
        disks: "local-disk " + select_first([DiskGB, auto_disk_size]) + " SSD"
        preemptible: 3
        maxRetries: 2
    }

    output {
        File VariantBed = "~{Prefix}.individual_variants.bed"
    }
}
