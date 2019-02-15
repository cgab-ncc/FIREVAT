# FIREVAT Main Functions
#
# Last revised date:
#   February 15, 2019
#
# Authors:
#   Andy Jinseok Lee (jinseok.lee@ncc.re.kr)
#   Hyunbin Kim (khb7840@ncc.re.kr)
#   Bioinformatics Analysis Team, National Cancer Center of Korea


RunFIREVATfromConfig <- function(vcf.file, 
                                 config.file,
                                 vcf.file.genome, 
                                 df.ref.mut.sigs, 
                                 target.mut.sigs, 
                                 sequencing.artifact.mut.sigs,
                                 num.cores, 
                                 output.dir, 
                                 pop.size = 200,
                                 max.iter = 200,
                                 run = 50, 
                                 pmutation = 0.25,
                                 verbose = TRUE)  {
    # Filters point mutations in the specified vcf.file based on mutational signature 
    # decomposition and outputs the filtered vcf data and metadata related to the
    # filtering process
    #
    # Args:
    #   vcf.file = character vector (path to the vcf.file)
    #   vcf.file.genome = reference genome of the vcf file (hg19 or hg38)
    #   vcf.filter = list (firevat_vcf::MakeMutect2Filter)
    #   vcf.calling.method = string (one of constants::Mutation.Calling.Methods)
    #   df.ref.mut.sigs = data.frame of reference mutational signatures
    #   sequencing.artifact.mut.sigs = character vector 
    #                                  (signature names used to perform filtering)
    #   optim.method = string ("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", or "Brent")
    #   num.cores = integer (number of cores available for the optimization computation)
    #   output.dir = string (output directory to which all resulting reports will be saved)
    #
    # Returns:
    #

    # 
    start.datetime <- Sys.time()
    if (verbose) {
        print(start.datetime)
        print("Initializing FIREVAT variant filtering pipeline")
    }

    # Create the output directory
    if (!dir.exists(output.dir))  {
        dir.create(output.dir, recursive = T)
    }
    else  {
        print("output.dir already exists")
    }
    
    # Read the vcf file
    vcf.obj <- ReadVCF(vcf.file = vcf.file,
                       genome = vcf.file.genome)

    config.obj <- ParseConfigFile(config.file, verbose = verbose)
    
    # Initially select point mutations with all 7 FORMAT parameters and consider
    # this the initial set of mutations
    vcf.objs <- InitializeVCFFromConfig(vcf.obj = vcf.obj,
                                        config.obj = config.obj,
                                        verbose = verbose)

    vcf.obj <- vcf.objs$vcf.obj.filtered


    vcf.filter <- MakeFilterFromConfig(config.obj)

    if (nrow(vcf.obj$data) <= 50)  {
        print("FIREVAT must have at least 50 mutations to run.")
        return()
    }
    
    print(paste0("Starting with ", nrow(vcf.obj$data), " point mutations"))

    bits.list <- ParameterToBits(vcf.obj, config.obj, vcf.filter) 
    params.bit.len <- bits.list$params.bit.len
    vcf.obj <- bits.list$vcf.obj

    n.bits <- sum(params.bit.len)
    
    # Prepare data for optimization
    data = list(n.bits = n.bits,
                params.bit.len = params.bit.len,
                vcf.file = vcf.file,
                vcf.file.basename = gsub("\\.vcf", "", basename(vcf.file)),
                vcf.obj = vcf.obj,
                config.file = config.file,
                config.obj = config.obj,
                vcf.filter = vcf.filter,
                df.ref.mut.sigs = df.ref.mut.sigs,
                df.ref.mut.sigs.cosmic = GetCOSMICMutSigs(),
                target.mut.sigs = target.mut.sigs,
                sequencing.artifact.mut.sigs = sequencing.artifact.mut.sigs,
                output.dir = output.dir,
                start.datetime = start.datetime,
                pmutation = pmutation,
                pop.size = pop.size,
                max.iter = max.iter,
                run = run)

    # Optimize filter parameters
    print("Started running optimization part")
    ga.results <- ga(type = "binary",
                     fitness =  function(string) GAOptimizationObjFn(string, data), 
                     nBits = data$n.bits,  
                     popSize = pop.size, 
                     maxiter = max.iter,
                     parallel = num.cores,
                     run = run,
                     pmutation = pmutation,
                     monitor = function(obj) GAMonitorFn(obj, data),
                     keepBest = T)
    # parallel = num.cores,
    
    print(summary(ga.results))
    print("Finished running optimization part")
    
    x.solution.binary <- as.numeric(ga.results@solution[1,]) 
    
    x.solution.decimal <- GADecodeBinaryString(
        string = x.solution.binary, data = data
    )
    
    print("Optimized binary solution")
    print(x.solution.binary)
    
    print("Optimized decimal solution")
    print(x.solution.decimal)
    
    data$end.datetime <- Sys.time()
    
    # Report results
    ReportFIREVATResultsfromConfig(x.solution.decimal = x.solution.decimal,
                                   data = data)

    return(list(x.solution.decimal = x.solution.decimal, data = data))
}
