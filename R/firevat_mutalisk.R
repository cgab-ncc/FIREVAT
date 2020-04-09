# FIREVAT Mutalisk Functions
# Adapted from Mutalisk (Lee et al., Nucleic Acids Research 2018; PMID 29790943)
#
# Last revised date:
#   February 19, 2019
#
# Authors:
#   Andy Jinseok Lee (jinseok.lee@ncc.re.kr)
#   Hyunbin Kim (khb7840@ncc.re.kr)
#   Bioinformatics Analysis Team, National Cancer Center Korea


#' @title MutaliskParseVCFObj
#' @description Parses a vcf.obj and prepares it to run Mutalisk.
#'
#' @param vcf.obj A list from ReadVCF
#'
#' @return A data.frame
#'
#' @importFrom deconstructSigs mut.to.sigs.input
#' @importFrom BSgenome getBSgenome
#' @export
MutaliskParseVCFObj <- function(vcf.obj)  {

    bsg <- BSgenome::getBSgenome(vcf.obj$genome)

    vcf.obj$data$Sample <- rep("sample", nrow(vcf.obj$data))
    df.deconstructsigs.sigs.input <- mut.to.sigs.input(mut.ref = vcf.obj$data,
                                                       sample.id = "Sample",
                                                       chr = "CHROM",
                                                       pos = "POS",
                                                       ref = "REF",
                                                       alt = "ALT",
                                                       bsg = bsg)

    df.temp <- t(df.deconstructsigs.sigs.input)

    mutation.types <- rownames(df.temp)
    context.both.strand <- c()
    for (curr.mutation.type in mutation.types)  {
        ref.1 <- sub('(.*)\\[.*>.*\\].*', '\\1', curr.mutation.type)
        alt.2 <- sub('.*\\[(.*)>.*\\].*', '\\1', curr.mutation.type)
        alt.3 <- sub('.*\\[.*>(.*)\\].*', '\\1', curr.mutation.type)
        ref.4 <- sub('.*\\[.*>.*\\](.*)', '\\1', curr.mutation.type)
        context.both.strand <- c(context.both.strand, paste0(ref.1, alt.2, ">", alt.3, ref.4))
    }
    context.both.strand <- c(context.both.strand, "TOTAL")
    freq <- c(as.numeric(df.temp), sum(as.numeric(df.temp)))

    df.mutalisk.input <- data.frame(list(context_both_strand = context.both.strand,
                                         freq = freq),
                                    stringsAsFactors = F)
    return(df.mutalisk.input)
}


#' @title MutaliskSigNamesToIndices
#' @description Converts target.mut.sigs to an integer vector corresponding to the indices in df.ref.mut.sigs
#'
#' @param df.ref.mut.sigs A data.frame of reference mutational signatures
#' @param target.mut.sigs A character vector of target mutational signatures names
#'
#' @return An integer vector of target mutational signatures indices
#' (the first target mutational signature (e.g. 'SBS1') gets assigned a value of 1)
#'
#' @keywords internal
MutaliskSigNamesToIndices <- function(df.ref.mut.sigs,
                                      target.mut.sigs)  {
    target.mut.sigs.indices <- which(colnames(df.ref.mut.sigs) %in% target.mut.sigs)
    return(target.mut.sigs.indices - 3)
}


#' @title MutaliskSigIndicesToNames
#' @description Converts target mutational signatures indices to the corresponding target
#' mutational signatures integer vector corresponding to the indices in df.ref.mut.sigs
#'
#' @param df.ref.mut.sigs A data.frame of reference mutational signatures
#' @param target.mut.sigs.indices An integer vector of target mutational signatures indices
#' (the first PCAWG signature, SBS1, has a value of 1)
#'
#' @return A character vector of target mutational signatures names
#'
#' @keywords internal
MutaliskSigIndicesToNames <- function(df.ref.mut.sigs,
                                      target.mut.sigs.indices)  {
    target.mut.sigs.indices <- target.mut.sigs.indices + 3
    return(colnames(df.ref.mut.sigs)[target.mut.sigs.indices])
}


#' @title RunMutaliskHelper
#' @description Helper function for RunMutalisk
#'
#' @param vcf.trinucleotide.data A data.frame (from firevat_mutalisk::MutaliskParseVCFObj)
#' @param df.ref.mut.sigs A data.frame of reference mutational signatures
#' @param target.mut.sigs A character vector of target mutational signatures names
#'
#' @return A list with the following elements
#' \itemize{
#'  \item{num.point.mutations}{An integer value - count of total point mutations}
#'  \item{sub.types}{A character vector of length 96}
#'  \item{sub.types.spectrum}{A numeric vector of length 96}
#'  \item{num.mut.sigs}{An integer value (count of unique mutational signatures identified)}
#'  \item{identified.mut.sigs}{A character vector where each element is a mutational signature identified}
#'  \item{identified.mut.sigs.probs}{A numeric vector where each element is the weight of mutational signature identified.
#'   The ordering follows identified.mut.sigs}
#'  \item{identified.mut.sigs.spectrum}{A numeric vector of length 96}
#'  \item{residuals}{A numeric vector of length 96}
#'  \item{rss}{A numeric value (residual sum of squares)}
#'  \item{cos.sim.score}{A numeric value (cosine similarity score between observed mutational spectrum and
#'   reconstructed mutational signatures)}
#'  \item{all.models.sigs}{A list where each element is a model; a model is a list of signatures identified)}
#'  \item{all.models.sigs.probs}{A list where each element is a model; a model is a list of contribution probabilities}
#'  \item{all.models.cos.sim.scores}{A list where each element is a model; a model is a list of cosine similarity socres}
#' }
#'
#' @export
#' @import lsa
#' @importFrom caTools combs
RunMutaliskHelper <- function(vcf.trinucleotide.data,
                              df.ref.mut.sigs,
                              target.mut.sigs) {
    target_signature <- MutaliskSigNamesToIndices(df.ref.mut.sigs = df.ref.mut.sigs,
                                                  target.mut.sigs = target.mut.sigs)
    num_rounds <- length(target_signature)
    if(num_rounds > 7) { # get up to 7 models
        num_rounds <- 7
    }

    min_prop = 0.01
    zeta_value = 1e-10

    spec1 <- as.matrix(vcf.trinucleotide.data)
    sigdata <- df.ref.mut.sigs

    sigdata <- sigdata[order(sigdata$`Type`),]
    data1 <- cbind(spec = as.integer(spec1[1:96,2]),
                   sigdata[1:96,4:ncol(sigdata)])

    num_mut = sum(data1[,1])
    this_subclass = spec1[1:96,1]
    this_spectrum = data1$spec/num_mut

    zero <- which(data1[,1] == 0)
    data1 <- data1[data1[,1] != 0,]

    total_num <- as.integer(spec1[97,2])
    total_CA <- sum(as.integer(spec1[1:16,2] ))
    total_CG <- sum(as.integer(spec1[17:32,2] ))
    total_CT <- sum(as.integer(spec1[33:48,2] ))
    total_TA <- sum(as.integer(spec1[49:64,2] ))
    total_TC <- sum(as.integer(spec1[65:80,2] ))
    total_TG <- sum(as.integer(spec1[81:96,2] ))
    prop_CA = round(total_CA/total_num, 3)
    prop_CG = round(total_CG/total_num, 3)
    prop_CT = round(total_CT/total_num, 3)
    prop_TA = round(total_TA/total_num, 3)
    prop_TC = round(total_TC/total_num, 3)
    prop_TG = round(total_TG/total_num, 3)

    LIK <- c()
    bdb <- c()
    sdb <- list()
    pdb <- list()
    min_sig = c()
    min_par = c()
    min_result = c()
    stopIdx <- 0
    stopVal <- 0
    farr = c()

    for (nrs in 1:num_rounds) {
        min_value = 1000000000000000

        if (nrs == 1) {
            all_comb <- combs(target_signature,1)
            length_all <- length(all_comb[,1])
            for (this_n in 1:length_all) {
                background <- data1[,all_comb[this_n,]+1]+zeta_value
                ## BIC: Bayesian Information Criterion
                this_value <- sum(log(background)*data1[,1])*-2+log(num_mut)*1
                if (this_value < min_value) {
                    min_sig = all_comb[this_n,]
                    min_value = this_value; min_par = 1
                }
                LIK[nrs] <- min_value
            }
            stopVal <- LIK[nrs]
        } else {
            str <- "function(par, data, sigs) { sum((data[,1]-data[,sigs[1]]*par[1]"
            for (idx in 2:nrs) {
                if (idx == nrs) {
                    break
                }
                str <- paste0(str," - data[,sigs[",idx,"]]*par[",idx,"]")
            }
            str2 <- "par[1]"
            for (idx in 1:(nrs-1)) {
                if ((idx+1) == nrs) {
                    break
                }
                str2 <- paste0(str2,"-par[",(idx+1),"]")
            }
            str <- paste0(str," - data[,sigs[",nrs,"]]*max(0,num_mut-",str2,"))^2) }")
            min.RSS <- eval(parse(text = str)) ## Linear Regression

            # FAST VERSION
            if (nrs >= 4) {
                tmp = min_sig
                tmp <- tmp[which(!tmp %in% farr)]
                tmp_par <- min_par[which(!min_sig %in% farr)]
                max.idx <- which.max(tmp_par)
                fixed <- tmp[max.idx]
                farr <- c(farr, fixed)

                if (length(target_signature) > 40) {
                    if (nrs != length(target_signature)) {
                        new_target <- target_signature[!target_signature %in% farr]

                        ac <- combs(new_target, 3)
                        for (i in 1:length(farr)) {
                            ac <- as.matrix(cbind(farr[i],ac))
                        }
                        all_comb <- ac
                        length_all <- length(all_comb[,1])
                    } else {
                        all_comb <- combs(target_signature, nrs)
                        length_all <- length(all_comb[,1])
                    }
                } else {
                    if (length(farr) == 1) {
                        tmp = min_sig
                        tmp <- tmp[which(!tmp %in% farr)]
                        tmp_par <- min_par[which(!min_sig %in% farr)]
                        max.idx <- which.max(tmp_par)
                        fixed <- tmp[max.idx]
                        farr <- c(farr, fixed)
                    }

                    if (nrs != length(target_signature)) {
                        all_comb <- combs(target_signature, nrs)
                        us_combs <- all_comb
                        for (i in farr) {
                            str <- paste0("us_combs <- us_combs[apply(us_combs,1,function(x) any (x==",i,")),]")
                            eval(parse(text = str))
                        }
                        all_comb <- us_combs
                        length_all <- length(all_comb[,1])
                    } else {
                        all_comb <- combs(target_signature, nrs)
                        length_all <- length(all_comb[,1])
                    }
                }
            } else {
                all_comb <- combs(target_signature, nrs)
                length_all <- length(all_comb[,1])
            }

            for (this_n in 1:length_all) {
                result <- optim(par = rep(1, (nrs-1)),
                                fn = min.RSS,
                                data = data1,
                                sigs = c(all_comb[this_n,]) + 1,
                                method = "L-BFGS-B",
                                lower = min_prop*num_mut,
                                upper = num_mut)

                if (result$value < min_value) {
                    min_value = result$value
                    min_result = result
                    min_sig = c(all_comb[this_n,])
                }
            }
            min_par = c(min_result$par,max(0,num_mut-sum(min_result$par)))
            min_par = min_par/sum(min_par)

            sss <- ""
            for (i in 1:nrs) {
                sss <- paste0(sss,"data1[,min_sig[",i,"]+1]*min_par[",i,"]+")
            }
            sss <- paste0(sss,"zeta_value")
            background <- eval(parse(text = sss))

            ## BIC: Bayesian Information Criterion
            LIK[nrs] <- sum(log(background)*data1[,1])*-2+log(num_mut)*nrs

            if (nrs == 1) {
                stopVal <- LIK[nrs]
            } else {
                if (LIK[nrs] < stopVal) {
                    if (stopIdx == 1 ) {
                        stopIdx <- stopIdx - 1
                    }
                    stopVal <- LIK[nrs]
                } else {
                    stopIdx <- stopIdx + 1
                }
            }
        }

        sdb[nrs] <- list(min_sig)
        pdb[nrs] <- list(min_par)
        bdb <- cbind(bdb,background)

        if (stopIdx == 2) {
            if (nrs < 4) {
                stopIdx <- 0
            } else {
                num_rounds <- nrs
                break
            }
        }
    }

    for (idx in num_rounds:1) {
        if (min(LIK) == LIK[idx]) {
            final_signum = idx
            final_sig = sdb[[idx]]
            final_par = pdb[[idx]]
            final_background = bdb[,idx]
        }
    }

    if (length(final_background) < 96) {
        for(i in zero) {
            dtmp <- as.vector(final_background)
            if( i == 96 ) {
                str <- paste0("final_background <- c(dtmp[1:",i-1,"],0)")
            } else if (i == 1) {
                str <- paste0("final_background <- c(0,dtmp[",i,":",length(dtmp),"])")
            } else {
                str <- paste0("final_background <- c(dtmp[1:",i-1,"],0,dtmp[",i,":",length(dtmp),"])")
            }
            eval(parse(text = str))
        }
    }

    # Get cosine similarity for all models considered
    all.models.cos.sim.scores <- c()
    for (idx in 1:num_rounds) {
        curr.background = bdb[,idx]
        if (length(curr.background) < 96) {
            for(i in zero) {
                dtmp <- as.vector(curr.background)
                if (i == 96) {
                    str <- paste0("curr.background <- c(dtmp[1:",i-1,"],0)")
                } else if (i == 1) {
                    str <- paste0("curr.background <- c(0,dtmp[",i,":",length(dtmp),"])")
                } else {
                    str <- paste0("curr.background <- c(dtmp[1:",i-1,"],0,dtmp[",i,":",length(dtmp),"])")
                }
                eval(parse(text = str))
            }
        }

        curr.residuals = this_spectrum - curr.background
        curr.cos.similarity <- lsa::cosine(this_spectrum, curr.background)
        all.models.cos.sim.scores <- c(all.models.cos.sim.scores, curr.cos.similarity)
    }

    # Get signature names
    all.models.signatures <- list()
    curr.model.idx <- 1
    for (curr.model in sdb) {
        all.models.signatures[[curr.model.idx]] <- MutaliskSigIndicesToNames(df.ref.mut.sigs = df.ref.mut.sigs,
                                                                             target.mut.sigs.indices = curr.model)
        curr.model.idx <- curr.model.idx + 1
    }

    final_residuals = this_spectrum - final_background
    cosine_similarity <- lsa::cosine(this_spectrum, final_background)

    rss <- final_residuals * final_residuals
    rss <- sum(rss)

    return(list(num.point.mutations = num_mut,
                sub.types = this_subclass,
                sub.types.spectrum = this_spectrum,
                num.mut.sigs = final_signum,
                identified.mut.sigs = MutaliskSigIndicesToNames(
                    df.ref.mut.sigs = df.ref.mut.sigs,
                    target.mut.sigs.indices = final_sig),
                identified.mut.sigs.probs = final_par,
                identified.mut.sigs.spectrum = final_background,
                residuals = final_residuals,
                rss = rss,
                cos.sim.score = cosine_similarity,
                all.models.sigs = all.models.signatures,
                all.models.sigs.probs = pdb,
                all.models.cos.sim.scores = all.models.cos.sim.scores))
}


#' @title RunMutalisk
#' @description Identifies mutational signatures using Mutalisk
#'
#' @param vcf.obj A list (from firevat_vcf::ReadVCF)
#' @param df.ref.mut.sigs A data.frame of reference mutational signatures
#' @param target.mut.sigs A character vector of target mutational signatures names to identify from
#' @param random.sampling.candidate.mut.sigs A character vector of mutational signatures names
#' that gets appended to the list of candidate mutational signatures so that these are
#' always considered.
#' @param method A string value (must be either 'random.sampling' or 'all').
#' The method 'random.sampling' samples (without replacement) 'n.sample' number of signatures
#' 'n.iter' number of times and runs the candidate signatures one last time.
#' The method 'all' uses all target.mut.sigs
#' @param n.sample An integer value ('random.sampling' method parameter)
#' Number of signatures to choose for each iteration of random sampling).
#' @param n.iter An integer value ('random.sampling' method parameter).
#' Number of iterations to perform random sampling.
#' @param verbose If true, provides process details
#'
#' @return A list with the following elements
#' \itemize{
#'  \item{num.point.mutations}{An integer value - count of total point mutations}
#'  \item{sub.types}{A character vector of length 96}
#'  \item{sub.types.spectrum}{A numeric vector of length 96}
#'  \item{num.mut.sigs}{An integer value (count of unique mutational signatures identified)}
#'  \item{identified.mut.sigs}{A character vector where each element is a mutational signature identified}
#'  \item{identified.mut.sigs.probs}{A numeric vector where each element is the weight of mutational signature identified.
#'   The ordering follows identified.mut.sigs}
#'  \item{identified.mut.sigs.spectrum}{A numeric vector of length 96}
#'  \item{residuals}{A numeric vector of length 96}
#'  \item{rss}{A numeric value (residual sum of squares)}
#'  \item{cos.sim.score}{A numeric value (cosine similarity score between observed mutational spectrum and
#'   reconstructed mutational signatures)}
#'  \item{all.models.sigs}{A list where each element is a model; a model is a list of signatures identified)}
#'  \item{all.models.sigs.probs}{A list where each element is a model; a model is a list of contribution probabilities}
#'  \item{all.models.cos.sim.scores}{A list where each element is a model; a model is a list of cosine similarity socres}
#' }
#'
#' @export
#' @import lsa
RunMutalisk <- function(vcf.obj,
                        df.ref.mut.sigs,
                        target.mut.sigs,
                        random.sampling.candidate.mut.sigs = c(),
                        method = "random.sampling",
                        n.sample = 20,
                        n.iter = 10,
                        verbose = TRUE) {
    if (method != "all" && method != "random.sampling") {
        stop("The parameter 'method' must be either 'all' or 'random.sampling'")
    }

    # Parse vcf.obj and get trinucleotide data
    vcf.trinucleotide.data <- MutaliskParseVCFObj(vcf.obj = vcf.obj)

    # Identify mutational signatures among the target signatures
    if (verbose == TRUE) {
        PrintLog("* Started running Mutalisk with:")
        PrintLog(paste0("** ", length(target.mut.sigs), " target signatures"))
        PrintLog(paste0("** '", method, "' method"))
    }

    # Run Mutalisk based on 'method'
    if (method == "random.sampling")  {
        target.mut.sigs.sampled <- c()
        for (i in 1:n.iter)  {
            if (verbose == TRUE) {
                PrintLog(paste0("* Random sampling iteration: ", i))
            }
            candidate.mut.sigs <- sample(target.mut.sigs, n.sample)
            candidateresults <- RunMutaliskHelper(vcf.trinucleotide.data = vcf.trinucleotide.data,
                                                  df.ref.mut.sigs = df.ref.mut.sigs,
                                                  target.mut.sigs = candidate.mut.sigs)
            target.mut.sigs.sampled <- c(target.mut.sigs.sampled,
                                         candidateresults$identified.mut.sigs)
        }
        target.mut.sigs.sampled <- unique(c(target.mut.sigs.sampled, random.sampling.candidate.mut.sigs))
        if (verbose == TRUE) {
            PrintLog(paste0("* Running Mutalisk with ", length(target.mut.sigs.sampled), " candidate signatures"))
        }
        results <- RunMutaliskHelper(vcf.trinucleotide.data = vcf.trinucleotide.data,
                                     df.ref.mut.sigs = df.ref.mut.sigs,
                                     target.mut.sigs = target.mut.sigs.sampled)
    }
    if (method == "all")  {
        results <- RunMutaliskHelper(vcf.trinucleotide.data = vcf.trinucleotide.data,
                                     df.ref.mut.sigs = df.ref.mut.sigs,
                                     target.mut.sigs = target.mut.sigs)
    }

    if (verbose == TRUE) {
        PrintLog("* Finished running Mutalisk")
    }
    return(results)
}
