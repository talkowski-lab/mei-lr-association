version 1.0

import "utils/query_vcf.wdl" as QueryVCF
import "utils/gzip_file.wdl" as GzipFile

workflow FilterMEIVCF {
    input {
        File VCF
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

    call GzipFile.Gzip as Gzip {
      input:
          InputFile = VCFQuery.Output
    }

    output {
        File Output = Gzip.Output
    }
}
