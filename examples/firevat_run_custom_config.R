library(FIREVAT)


# Prepare files
sample.vcf.file <- system.file("extdata", "DCC_PCAWG_Cell_Lines_HCC1143.vcf", package = "FIREVAT")
config.file <- system.file("config", "PCAWG_DKFZ_Cell_Line_Filtering_Params.json", package = "FIREVAT")


# Run FIREVAT
results <- RunFIREVAT(vcf.file = sample.vcf.file,
                      vcf.file.genome = 'hg19',
                      config.file = config.file,
                      df.ref.mut.sigs = GetPCAWGMutSigs(),
                      target.mut.sigs = GetPCAWGMutSigsNames(),
                      sequencing.artifact.mut.sigs = PCAWG.All.Sequencing.Artifact.Signatures,
                      output.dir = "/home/jinseoklee/Documents/Projects/FIREVAT_Workspace/20190219_090000/",
                      num.cores = 6,
                      pop.size = 200,
                      max.iter = 10,
                      run = 50,
                      pmutation = 0.25)
