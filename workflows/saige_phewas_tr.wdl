version 1.0

workflow SAIGE_phewas_TR {
    meta {
        description: "SAIGE PheWAS of either tandem-repeat lengths (mode tr: FORMAT/DS scaled to [0, 2]) or ordinary biallelic variants (mode variant: 0/1/2 GT, or DS); tests one phecode per run."
    }

    parameter_meta {
        plink_bed_file: "PLINK genotypes used only to fit the null model and estimate variance ratios (Step 1); typically LD-pruned common SNPs. Sample IDs must match person_id. Required unless both precomputed files are given."
        precomputed_null_model_rda: "Optional .rda from an earlier Step 1 of the SAME phecode, covariates, samples and sparse GRM. Alone: Step 1 skips the model fit and only estimates variance ratios (--skipModelFitting). With precomputed_variance_ratio_txt: Step 1 is skipped entirely."
        precomputed_variance_ratio_txt: "Optional .varianceRatio.txt from the same earlier Step 1 as precomputed_null_model_rda; only valid together with it."
        sparseGRM: "Sparse GRM .mtx from SAIGE createSparseGRM.R; used in both steps."
        sparseGRM_IDlist: "Sample-ID file written alongside sparseGRM."
        phenotype_file: "Whitespace-delimited table with a person_id column, one column per phecode, and every covariate column."
        phecode: "The one phenotype column of phenotype_file to test in this run."
        covariate_list_file: "Text file, one covariate column name per line. Omit to use PC1-PC16, sex, age, age2, age_sex, age2_sex."
        categorical_covariate_list_file: "Text file, one categorical covariate per line; each must also be in covariate_list_file. Omit to use sex; pass an empty file for none."
        trait_type: "binary or quantitative; quantitative traits are inverse-normalized."
        mode: "tr (scaled repeat lengths as dosages, no AF filter, mean imputation) or variant (ordinary variants, standard AF/MAC filters, best-guess imputation). Selects which tr_* / variant_* filter set is used."
        vcf: "bgzipped VCF, one biallelic record per variant. mode tr: FORMAT/DS (Number=1, Float) = repeat length scaled to [0, 2], ALT can be a placeholder such as <TR>; effect sizes are per unit of scaled length. mode variant: hard-called FORMAT/GT, or FORMAT/DS dosages if vcffield = DS."
        vcf_csi: "CSI index of vcf (bcftools index -c)."
        vcffield: "FORMAT field SAIGE reads. Defaults to DS in mode tr and GT in mode variant."
        chromosome_list_file: "Text file, one chromosome per line, named exactly as in vcf (chr1 vs 1); each must have at least one record."
        ImageTag: "Tag of ayenkin1871/mei-lr-association-saige (latest or a 7-char commit SHA)."
        preflight_records: "Number of leading VCF records whose values (dosage range or GT alleles) are checked before Step 1."
    }

    input {
        File? plink_bed_file
        File? plink_bim_file
        File? plink_fam_file
        File? precomputed_null_model_rda
        File? precomputed_variance_ratio_txt
        File sparseGRM
        File sparseGRM_IDlist
        String phecode
        File vcf
        File vcf_csi
        File chromosome_list_file
        File? covariate_list_file
        File? categorical_covariate_list_file
        File phenotype_file
        String trait_type
        String output_prefix
        String ImageTag = "latest"
        String mode = "tr"
        String? vcffield
        Int preflight_records = 2000

        String memory = "8G"

        # SAIGE derives AF as mean(DS)/2, which for a scaled repeat length measures where the
        # population sits in the length range, not how rare an allele is. A TR whose lengths
        # cluster near either end of the scale would be dropped by a MAF filter despite real
        # variation, so the AF filter is off by default.
        Float tr_minimal_af = 0
        Float tr_min_mac = 0.5
        Float tr_max_missing = 0.5

        Float variant_minimal_af = 0.0001
        Float variant_min_mac = 20
        Float variant_max_missing = 0.15

        # MAC category bounds for categorical variance ratios, shared by Step 1 and Step 2 so the
        # two cannot drift. SAIGE requires length(min_exclude) == length(max_include) + 1, and
        # length(min_exclude) == the number of "null" rows Step 1 writes to the varianceRatio file.
        # If these disagree between steps, SAIGE applies each ratio to the wrong MAC bin silently.
        String cate_var_ratio_min_mac_exclude = "100,500"
        String cate_var_ratio_max_mac_include = "500"

        Int cpu = 2
        # Empty Step 1 outputs are often a full /tmp or a full local-disk, not a SAIGE bug.
        String disk_size = "local-disk 200 HDD"
        Int preemptible = 1
    }

    String saige_docker = "ayenkin1871/mei-lr-association-saige:" + ImageTag

    Boolean is_tr = mode == "tr"
    String field = select_first([vcffield, if is_tr then "DS" else "GT"])
    Float minimal_af = if is_tr then tr_minimal_af else variant_minimal_af
    Float min_mac = if is_tr then tr_min_mac else variant_min_mac
    Float max_missing = if is_tr then tr_max_missing else variant_max_missing
    # SAIGE's default best_guess rounds a missing value to round(2*AF), right for hard calls
    # but not for a continuous length.
    String impute_method = if is_tr then "mean" else "best_guess"

    Boolean skip_step1 = defined(precomputed_null_model_rda) && defined(precomputed_variance_ratio_txt)

    # read_lines keeps blank lines and trailing \r; the tasks drop both when they word-split.
    Array[String] chromosome_list = read_lines(chromosome_list_file)
    Array[String] covariate_list = if defined(covariate_list_file)
        then read_lines(select_first([covariate_list_file]))
        else ["PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","sex","age","age2","age_sex","age2_sex"]
    Array[String] categorical_covariate_list = if defined(categorical_covariate_list_file)
        then read_lines(select_first([categorical_covariate_list_file]))
        else ["sex"]

    call ValidateVCF {
        input:
            vcf               = vcf,
            vcf_csi           = vcf_csi,
            mode              = mode,
            vcffield          = field,
            chromo_list       = chromosome_list,
            phenotype_file    = phenotype_file,
            phecode           = phecode,
            covariate_list    = covariate_list,
            categorical_covariate = categorical_covariate_list,
            n_records         = preflight_records,
            has_plink         = defined(plink_bed_file) && defined(plink_bim_file) && defined(plink_fam_file),
            has_null_model    = defined(precomputed_null_model_rda),
            has_variance_ratio = defined(precomputed_variance_ratio_txt),
            saige_docker      = saige_docker
    }

    if (!skip_step1) {
        # select_first only runs after ValidateVCF has confirmed the PLINK files were given.
        call RunFitNullGLMM {
            input:
                preflight_report  = ValidateVCF.report,
                existing_rda      = precomputed_null_model_rda,
                plink_bed_file    = select_first([plink_bed_file]),
                plink_bim_file    = select_first([plink_bim_file]),
                plink_fam_file    = select_first([plink_fam_file]),
                sparseGRM         = sparseGRM,
                sparseGRM_IDlist  = sparseGRM_IDlist,
                phenotype_file    = phenotype_file,
                phecode           = phecode,
                trait_type        = trait_type,
                output_prefix     = output_prefix,
                saige_docker      = saige_docker,
                cate_var_ratio_min_mac_exclude = cate_var_ratio_min_mac_exclude,
                cate_var_ratio_max_mac_include = cate_var_ratio_max_mac_include,
                covariate_list = covariate_list,
                categorical_covariate = categorical_covariate_list,

                memory            = memory,
                cpu               = cpu,
                disk              = disk_size,
                preemptible       = preemptible
        }
    }

    File model_rda = select_first([RunFitNullGLMM.null_model_rda, precomputed_null_model_rda])
    File variance_ratio = select_first([RunFitNullGLMM.variance_ratio_txt, precomputed_variance_ratio_txt])

    call RunStep2 {
        input:
            vcf              = vcf,
            vcf_csi          = vcf_csi,
            mode             = mode,
            impute_method    = impute_method,
            preflight_report = ValidateVCF.report,
            variance_ratio   = variance_ratio,
            model_file       = model_rda,
            sparseGRM        = sparseGRM,
            sparseGRM_IDlist = sparseGRM_IDlist,
            vcffield         = field,
            chromo_list      = chromosome_list,
            phecode          = phecode,
            output_prefix    = output_prefix,
            minimal_af       = minimal_af,
            min_mac          = min_mac,
            max_missing      = max_missing,
            cate_var_ratio_min_mac_exclude = cate_var_ratio_min_mac_exclude,
            cate_var_ratio_max_mac_include = cate_var_ratio_max_mac_include,
            memory           = memory,
            saige_docker     = saige_docker,
            cpu              = cpu,
            disk             = disk_size,
            preemptible      = preemptible
    }

    output {
        File        preflight_report    = ValidateVCF.report
        File        null_model_rda      = model_rda
        File        variance_ratio_txt  = variance_ratio
        Array[File] step2_outputs       = RunStep2.outputfiles
    }
}


# Fails in minutes on input problems that would otherwise surface only after hours of Step 1,
# or not at all (out-of-range dosages and absent samples do not make SAIGE error).
task ValidateVCF {

    input {
        File vcf
        File vcf_csi
        String mode
        String vcffield
        Array[String] chromo_list
        File phenotype_file
        String phecode
        Array[String] covariate_list
        Array[String] categorical_covariate
        Int n_records
        Boolean has_plink
        Boolean has_null_model
        Boolean has_variance_ratio
        String saige_docker
    }

    Int disk_gb = ceil(size(vcf, "GB") + size(phenotype_file, "GB")) + 20

    command <<<
        set -euo pipefail

        ln -s "~{vcf}" tr.vcf.gz
        ln -s "~{vcf_csi}" tr.vcf.gz.csi
        errors=0
        fail() { echo "ERROR: $*" | tee -a preflight.txt >&2; errors=$((errors + 1)); }
        note() { echo "$*" | tee -a preflight.txt >&2; }
        # Drops the blank entries and \r that read_lines keeps from list files.
        clean_list() { local x; for x in "$@"; do x="${x//$'\r'/}"; if [[ -n "${x}" ]]; then printf '%s\n' "${x}"; fi; done; }
        mapfile -t chromos < <(clean_list ~{sep=' ' chromo_list})
        mapfile -t covars < <(clean_list ~{sep=' ' covariate_list})
        mapfile -t qcovars < <(clean_list ~{sep=' ' categorical_covariate})
        [[ ${#chromos[@]} -gt 0 ]] || fail "chromosome list file is empty"

        case "~{mode}" in
            tr|variant) ;;
            *) fail "mode must be tr or variant, got '~{mode}'" ;;
        esac
        if ~{has_variance_ratio} && ! ~{has_null_model}; then
            fail "precomputed_variance_ratio_txt was given without precomputed_null_model_rda"
        fi
        if ! { ~{has_null_model} && ~{has_variance_ratio}; } && ! ~{has_plink}; then
            fail "PLINK bed/bim/fam are required unless both precomputed_null_model_rda and precomputed_variance_ratio_txt are given"
        fi
        if [[ "~{mode}" == "tr" && "~{vcffield}" == "GT" ]]; then
            fail "mode tr reads scaled lengths from a dosage field; GT is only valid in mode variant"
        fi

        hdr="$(bcftools view -h tr.vcf.gz | grep "^##FORMAT=<ID=~{vcffield}," || true)"
        if [[ -z "${hdr}" ]]; then
            fail "no FORMAT/~{vcffield} in the VCF header"
        elif [[ "${hdr}" != *"Number=1,"* ]]; then
            fail "FORMAT/~{vcffield} must be Number=1 (one value per sample); header: ${hdr}"
        fi

        multi="$( { bcftools view -H -m3 tr.vcf.gz || true; } | head -n 1 | cut -f1-5)"
        if [[ -n "${multi}" ]]; then
            fail "multi-allelic record found (${multi}); SAIGE needs one REF/ALT pair per TR"
        fi

        bcftools index -s tr.vcf.gz | cut -f1 | sort -u > vcf_contigs.txt
        for chromo in "${chromos[@]}"; do
            if ! grep -qxF "${chromo}" vcf_contigs.txt; then
                fail "chromosome ${chromo} has no records in the VCF (contigs present: $(paste -sd, vcf_contigs.txt))"
            fi
        done

        if [[ "~{vcffield}" == "GT" ]]; then
            { bcftools query -f '[%GT\t]\n' tr.vcf.gz || true; } | head -n ~{n_records} > gt_sample.txt
            read -r n_vals n_missing n_bad < <(awk -F'\t' '
                { for (i = 1; i <= NF; i++) {
                    if ($i == "") continue
                    n++
                    if ($i ~ /^\.([\/|]\.)?$/) { miss++; continue }
                    if ($i !~ /^[01.]([\/|][01.])?$/) bad++
                } }
                END { print n+0, miss+0, bad+0 }' gt_sample.txt)
            note "First ~{n_records} records: ${n_vals} GT calls, ${n_missing} missing"
            [[ "${n_vals}" -gt 0 ]] || fail "no GT values read from the VCF"
            [[ "${n_bad}" -eq 0 ]] || fail "${n_bad} GT calls are not 0/1 alleles (unsplit multi-allelic or malformed)"
        else
            { bcftools query -f '[%~{vcffield}\t]\n' tr.vcf.gz || true; } | head -n ~{n_records} > ds_sample.txt
            read -r n_vals n_missing n_bad n_out vmin vmax < <(awk -F'\t' '
                BEGIN { min = 1e300; max = -1e300 }
                { for (i = 1; i <= NF; i++) {
                    if ($i == "") continue
                    n++
                    if ($i == ".") { miss++; continue }
                    if ($i !~ /^-?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/) { bad++; continue }
                    v = $i + 0
                    if (v < 0 || v > 2) out++
                    if (v < min) min = v
                    if (v > max) max = v
                } }
                END { print n+0, miss+0, bad+0, out+0, (n > miss + bad ? min : "NA"), (n > miss + bad ? max : "NA") }' ds_sample.txt)
            note "First ~{n_records} records: ${n_vals} values, ${n_missing} missing, range [${vmin}, ${vmax}]"
            [[ "${n_vals}" -gt 0 ]] || fail "no ~{vcffield} values read from the VCF"
            [[ "${n_bad}" -eq 0 ]] || fail "${n_bad} non-numeric ~{vcffield} values"
            [[ "${n_out}" -eq 0 ]] || fail "${n_out} ~{vcffield} values outside [0, 2]; dosages must be in [0, 2]"
        fi

        bcftools query -l tr.vcf.gz | sort -u > vcf_samples.txt
        awk 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "person_id") c = i; if (!c) exit 1; next }
             { print $c }' "~{phenotype_file}" | sort -u > pheno_samples.txt \
            || fail "phenotype file has no person_id column"
        n_overlap="$(comm -12 vcf_samples.txt pheno_samples.txt | wc -l)"
        note "Samples: $(wc -l < vcf_samples.txt) in VCF, $(wc -l < pheno_samples.txt) in phenotype file, ${n_overlap} shared"
        [[ "${n_overlap}" -gt 0 ]] || fail "no VCF sample IDs match phenotype person_id"

        head -n 1 "~{phenotype_file}" | tr -d '\r' | tr -s ' \t' '\n' > pheno_columns.txt
        for covar in "${covars[@]}"; do
            grep -qxF "${covar}" pheno_columns.txt \
                || fail "covariate ${covar} is not a column of the phenotype file"
        done
        for qcovar in "${qcovars[@]}"; do
            printf '%s\n' "${covars[@]}" | grep -qxF "${qcovar}" \
                || fail "categorical covariate ${qcovar} is not in the covariate list"
        done
        grep -qxF "~{phecode}" pheno_columns.txt \
            || fail "phecode ~{phecode} is not a column of the phenotype file"
        # SAIGE pastes these names into an R formula, so tech-GC parses as tech minus GC.
        for name in "~{phecode}" "${covars[@]}"; do
            [[ "${name}" =~ ^([A-Za-z]|\.[A-Za-z._])[A-Za-z0-9._]*$ ]] \
                || fail "column name '${name}' is not a valid R name (letters, digits, . and _ only; must not start with a digit or _); rename it in the phenotype file"
        done

        if [[ "${errors}" -gt 0 ]]; then
            echo "${errors} preflight check(s) failed" >&2
            exit 1
        fi
        note "Preflight passed"
    >>>

    output {
        File report = "preflight.txt"
    }

    runtime {
        docker:      saige_docker
        memory:      "4G"
        cpu:         1
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: 1
    }
}


task RunFitNullGLMM {

    input {
        # Unused; makes Step 1 wait for ValidateVCF.
        File preflight_report
        File? existing_rda
        File plink_bed_file
        File plink_bim_file
        File plink_fam_file
        File sparseGRM
        File sparseGRM_IDlist
        File phenotype_file
        String phecode
        Array[String] covariate_list
        Array[String] categorical_covariate
        String trait_type
        String output_prefix
        String saige_docker

        # Must match the values passed to step2_SPAtests.R. See the workflow-level inputs.
        String cate_var_ratio_min_mac_exclude
        String cate_var_ratio_max_mac_include

        String memory
        Int cpu
        String disk
        Int preemptible
    }

    String inv_norm_arg = if (trait_type == "quantitative") then "--invNormalize=TRUE" else ""

    command <<<
        set -uxo pipefail

        # R tempfile()/save() default to /tmp, which on Terra/dataproc is a small
        # boot disk. A failed write leaves a 0-byte .rda behind.
        export TMPDIR="${PWD}/tmp" TMP="${PWD}/tmp" TEMP="${PWD}/tmp"
        mkdir -p "${TMPDIR}"
        df -h . /tmp >&2 || true

        # Drops the blank entries and \r that read_lines keeps from list files.
        clean_list() { local x; for x in "$@"; do x="${x//$'\r'/}"; if [[ -n "${x}" ]]; then printf '%s\n' "${x}"; fi; done; }
        covar_arg="$(clean_list ~{sep=' ' covariate_list} | paste -sd, -)"
        qcovar_arg="$(clean_list ~{sep=' ' categorical_covariate} | paste -sd, -)"
        model_file="~{output_prefix}_step1Out_~{phecode}.rda"
        vr_file="~{output_prefix}_step1Out_~{phecode}.varianceRatio.txt"

        # --skipModelFitting loads the model from exactly <outputPrefix>.rda.
        skip_fit_arg=""
        existing_rda="~{default='' existing_rda}"
        if [[ -n "${existing_rda}" ]]; then
            cp "${existing_rda}" "${model_file}"
            skip_fit_arg="--skipModelFitting=TRUE"
        fi

        # SAIGE often writes .rda + .varianceRatio.txt and then exits non-zero
        # (R warnings / glm.fit). Success is the files, not the R exit code.
        step1_fitNULLGLMM.R \
            --bedFile="~{plink_bed_file}" \
            --bimFile="~{plink_bim_file}" \
            --famFile="~{plink_fam_file}" \
            --useSparseGRMtoFitNULL=TRUE \
            --useSparseGRMforVarRatio=TRUE \
            --sparseGRMFile="~{sparseGRM}" \
            --sparseGRMSampleIDFile="~{sparseGRM_IDlist}" \
            --phenoFile="~{phenotype_file}" \
            --phenoCol="~{phecode}" \
            --covarColList="${covar_arg}" \
            --qCovarColList="${qcovar_arg}" \
            --sampleIDColinphenoFile=person_id \
            --traitType=~{trait_type} \
            --isCateVarianceRatio=TRUE \
            --cateVarRatioMinMACVecExclude=~{cate_var_ratio_min_mac_exclude} \
            --cateVarRatioMaxMACVecInclude=~{cate_var_ratio_max_mac_include} \
            --IsOverwriteVarianceRatioFile=TRUE \
            --outputPrefix="~{output_prefix}_step1Out_~{phecode}" \
            ${skip_fit_arg} \
            ~{inv_norm_arg} 2>&1 | tee step1.log
        rc=${PIPESTATUS[0]}
        # -s: exist AND non-empty. SAIGE creates 0-byte placeholders, then
        # save() can fail (full /tmp, full local-disk) and leave them empty.
        if [[ ! -s "${model_file}" || ! -s "${vr_file}" ]]; then
            echo "ERROR: step1 for ~{phecode} failed (exit ${rc}); ${model_file} $(stat -c%s "${model_file}" 2>/dev/null || echo 0) bytes, ${vr_file} $(stat -c%s "${vr_file}" 2>/dev/null || echo 0) bytes" >&2
            echo "---- SAIGE error lines ----" >&2
            grep -n -i -E 'error|stop|cannot|killed|no space' step1.log | tail -n 15 >&2 || true
            echo "---- last 40 lines of SAIGE output ----" >&2
            tail -n 40 step1.log >&2
            exit 1
        fi
        if [[ "${rc}" -ne 0 ]]; then
            echo "WARNING: step1 exited ${rc} but outputs exist; keeping" >&2
        fi
    >>>

    output {
        File null_model_rda     = "~{output_prefix}_step1Out_~{phecode}.rda"
        File variance_ratio_txt = "~{output_prefix}_step1Out_~{phecode}.varianceRatio.txt"
    }

    runtime {
        docker:      saige_docker
        memory:      memory
        cpu:         cpu
        disks:       disk
        preemptible: preemptible
    }
}


task RunStep2 {

    input {
        # Unused; makes Step 2 wait for ValidateVCF when Step 1 is skipped.
        File preflight_report
        File vcf
        File vcf_csi
        String mode
        String impute_method
        File model_file
        File variance_ratio
        # Step 1 always fits with --useSparseGRMtoFitNULL=TRUE, so Step 2 must be given the same
        # sparse GRM — otherwise SAIGE silently falls back to the "null" variance ratios.
        File sparseGRM
        File sparseGRM_IDlist
        String vcffield
        Array[String] chromo_list
        String phecode
        String output_prefix
        Boolean is_fastTest = true
        Float minimal_af
        Float min_mac
        Float max_missing
        String cate_var_ratio_min_mac_exclude
        String cate_var_ratio_max_mac_include

        String memory
        String saige_docker
        Int cpu
        String disk
        Int preemptible
    }

    command <<<
        set -euxo pipefail

        export TMPDIR="${PWD}/tmp" TMP="${PWD}/tmp" TEMP="${PWD}/tmp"
        mkdir -p "${TMPDIR}"

        # Drops the blank entries and \r that read_lines keeps from list files.
        clean_list() { local x; for x in "$@"; do x="${x//$'\r'/}"; if [[ -n "${x}" ]]; then printf '%s\n' "${x}"; fi; done; }
        mapfile -t chromo_list < <(clean_list ~{sep=' ' chromo_list})

        for chromo in "${chromo_list[@]}"; do
            step2_SPAtests.R \
                --vcfFile=~{vcf} \
                --vcfFileIndex=~{vcf_csi} \
                --vcfField=~{vcffield} \
                --AlleleOrder=ref-first \
                --chrom=${chromo} \
                --SAIGEOutputFile=~{output_prefix}_step2Out_~{phecode}_${chromo}_~{mode} \
                --minMAF=~{minimal_af} \
                --minMAC=~{min_mac} \
                --maxMissing=~{max_missing} \
                --impute_method=~{impute_method} \
                --sparseGRMFile="~{sparseGRM}" \
                --sparseGRMSampleIDFile="~{sparseGRM_IDlist}" \
                --cateVarRatioMinMACVecExclude=~{cate_var_ratio_min_mac_exclude} \
                --cateVarRatioMaxMACVecInclude=~{cate_var_ratio_max_mac_include} \
                --GMMATmodelFile="~{model_file}" \
                --varianceRatioFile="~{variance_ratio}" \
                --is_Firth_beta=TRUE \
                ~{true="--is_fastTest=TRUE" false="--is_fastTest=FALSE" is_fastTest} \
                --LOCO=FALSE \
                --is_output_moreDetails=TRUE
        done
    >>>

    output {
        Array[File] outputfiles = glob("~{output_prefix}_step2Out_~{phecode}_*_~{mode}*")
    }

    runtime {
        docker:      saige_docker
        memory:      memory
        cpu:         cpu
        disks:       disk
        preemptible: preemptible
    }
}
