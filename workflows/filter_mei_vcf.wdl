version 1.0

import "utils/query_vcf.wdl" as QueryVCF
import "utils/gzip_file.wdl" as GzipFile

workflow FilterMEIVCF {
    input {
        File VCF
        File MEIClusterBed
        String MEIType
    }

    String OutputName = MEIType + "_gt_mat.txt"
    String IncludeExpr = "INFO/SVTYPE=\"" + MEIType + "\""
    String QueryFormat = "%CHROM\t%POS\t%REF\t%ALT[\t%GT]"

    call QueryVCF.BcftoolsQuery as VCFQuery {
        input:
            VCF = VCF,
            QueryFormat = QueryFormat,
            IncludeExpr = IncludeExpr,
            OutputName = OutputName
    }

    call ClusterMEIGT {
        input:
            MEIGTMat = VCFQuery.Output,
            MEIClusterBed = MEIClusterBed,
            MEIType = MEIType
        
    }

    call GzipFile.Gzip as GzipGTMat {
      input:
          InputFile = VCFQuery.Output
    }

    call GzipFile.Gzip as GzipGTClustered {
      input:
          InputFile = ClusterMEIGT.ClusteredGT
    }

    output {
        File GTMat = GzipGTMat.Output
        File GTClustered = GzipGTClustered.Output
        File ExclClusters = ClusterMEIGT.ExclClusters
    }
}


task ClusterMEIGT {
    input {
      File MEIGTMat
      File MEIClusterBed
      String MEIType
      String ImageTag = "latest"
      Int MemoryGB = 4
      Int? DiskGB
    }

    Int auto_disk_size = ceil(size([MEIGTMat, MEIClusterBed], "GB") * 2) + 10

    String ClusteredGT_out = MEIType + "_gt_tidy_clustered.txt"
    String ExclClusters_out = MEIType + "_gt_tidy_dupclusters.txt"

    command <<<
        set -euo pipefail

        Rscript /scripts/cluster_mei_gt.R \
          --mei-gt-mat ~{MEIGTMat} \
          --mei-cluster-bed ~{MEIClusterBed} \
          --prefix ~{MEIType}
    >>>

    runtime {
        docker: "ayenkin1871/mei-lr-association-r_analysis:" + ImageTag
        memory: MemoryGB + " GB"
        cpu: 2
        disks: "local-disk " + select_first([DiskGB, auto_disk_size]) + " SSD"
        preemptible: 3
        maxRetries: 2
    }

    output {
      File ClusteredGT = ClusteredGT_out
      File ExclClusters = ExclClusters_out
    }

}

