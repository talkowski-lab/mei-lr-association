version 1.0

workflow ExtractVariants {
    input {
        String VDSUri
        File VariantList
        File SampleList
        String Branch = "main"
    }
    
    call ExtractVariantsVDS {
        input:
            VDSUri = VDSUri,
            VariantList = VariantList,
            SampleList = SampleList,
            Branch = Branch
    }

    output {
        File Genotypes = ExtractVariantsVDS.GenotypeTSV
    }
}

task ExtractVariantsVDS {
    input {
        String VDSUri
        File VariantList
        File SampleList
        String Branch = "main"
    }

    command <<<
        export SPARK_LOCAL_DIRS=/cromwell_root

        # writes VCF to bucket path 
        # and also generates outpath.txt upon completion 
        # of writing VCF 
        python3 /extract_snps.py \
            --variant-list ~{VariantList} \
            --sample-list ~{SampleList} \
            --gs-vds-uri ~{VDSUri} 
    >>>

    runtime {
        docker: "ghcr.io/alyenkin/mei-lr-association:" + Branch
        memory: "128G"
        cpu: 64
        disks: "local-disk 100 SSD"
    }
    
    output {
        File GenotypeTSV = "genotypes.tsv" 
    }
}
