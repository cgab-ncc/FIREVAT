library(FIREVAT)


# If you wish to run the sample code, change the 'output.dir' parameter
# and run the whole script.

# With 2 cores (each with 3.5GHz clock speed) this sample script
# takes about 10 minutes to complete

# Output directory
output.dir <- "pwd/home/jinseoklee/Documents/Projects/Temp/FIREVAT_TEMP_2/"

# Sample VCF file
sample.vcf.file <- system.file("extdata", "DCC_PCAWG_Cell_Lines_HCC1954.vcf", package = "FIREVAT")

# Config file
config.file <- system.file("config", "PCAWG_DKFZ_Cell_Line_Filtering_Params.json", package = "FIREVAT")

# Annotation DB
clinvar.vcf.file <- system.file("extdata", "clinvar_hg19_20190212.vcf", package = "FIREVAT")
clinvar.vcf.obj <- ReadVCF(vcf.file = clinvar.vcf.file, genome = "hg19", split.info = TRUE)
df.annotation.db <- PrepareAnnotationDB(annotation.vcf.obj = clinvar.vcf.obj)

# Annotation parameters
cols.to.display = c("GENEINFO", "CLNSIG")
filter.key.value.pairs <- list("CLNSIG" = c("Pathogenic", "Pathogenic/Likely_pathogenic", "Likely_pathogenic"))

# Run FIREVAT
results <- RunFIREVAT(vcf.file = sample.vcf.file,
                      vcf.file.genome = 'hg19',
                      config.file = config.file,
                      df.ref.mut.sigs = GetPCAWGMutSigs(),
                      target.mut.sigs = GetPCAWGMutSigsNames(),
                      sequencing.artifact.mut.sigs = PCAWG.All.Sequencing.Artifact.Signatures,
                      output.dir = output.dir,
                      objective.fn = Default.Obj.Fn,
                      num.cores = 4,
                      ga.pop.size = 200,
                      ga.max.iter = 100,
                      ga.run = 100,
                      ga.pmutation = 0.1,
                      perform.strand.bias.analysis = TRUE,
                      ref.forward.strand.var = "TumorDPRefForward",
                      ref.reverse.strand.var = "TumorDPRefReverse",
                      alt.forward.strand.var = "TumorDPAltForward",
                      alt.reverse.strand.var = "TumorDPAltReverse",
                      annotate = TRUE,
                      df.annotation.db = df.annotation.db,
                      annotated.columns.to.display = cols.to.display,
                      annotation.filter.key.value.pairs = filter.key.value.pairs,
                      annotation.filter.condition = "AND")
