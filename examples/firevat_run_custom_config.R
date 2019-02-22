library(FIREVAT)


# Sample VCF file
sample.vcf.file <- system.file("extdata", "DCC_PCAWG_Cell_Lines_HCC1954_GERMLINE.vcf", package = "FIREVAT")

# Config file
config.file <- system.file("config", "PCAWG_DKFZ_Cell_Line_Filtering_Params.json", package = "FIREVAT")

# Run FIREVAT
results <- RunFIREVAT(vcf.file = sample.vcf.file,
                      vcf.file.genome = 'hg19',
                      config.file = config.file,
                      df.ref.mut.sigs = GetPCAWGMutSigs(),
                      target.mut.sigs = GetPCAWGMutSigsNames(),
                      sequencing.artifact.mut.sigs = PCAWG.All.Sequencing.Artifact.Signatures,
                      output.dir = "",
                      num.cores = 6,
                      ga.pop.size = 200,
                      ga.max.iter = 100,
                      ga.run = 20,
                      ga.pmutation = 0.25,
                      mutalisk.random.sampling.count = 20,
                      mutalisk.random.sampling.max.iter = 10,
                      annotate = FALSE)

