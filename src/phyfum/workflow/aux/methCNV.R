if (!require("pacman")) install.packages("pacman")
p_load(optparse, cli)

# Define the command line options
option_list <- list(
  make_option(
    c("--input", "-i"),
    type = "character",
    default = NULL,
    help = "R object containing the methylation data, obtained with readIDAT.R (required)"
  ),
  make_option(
    c("--mc.cores", "-c"),
    type = "integer",
    default = 1,
    help = "Number of CPU cores to use [default: %default])"
  ),
  make_option(
    c("--genome_plot", "-g"),
    action="store_true", default=FALSE,
    help = "Print plots [default: %default]"
  )
)

# Parse the command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Libraries
cli_alert_info("Preparing the environment")

p_load(conumee, minfi, parallel, tibble, tidyr, dplyr, data.table, gtools)
p_load_gh("crukci-bioinformatics/rascal")

# Check if libraryd arguments are provided
if (is.null(opt$input)) {
  cat("Error: Please provide an input .Rdata.\n")
  cat("Run with --help for usage information.\n")
  quit(status = 1)
}

load(opt$input)
mc.cores <- opt$mc.cores
genome_plot <- opt$genome_plot
RNGkind("L'Ecuyer-CMRG") #Different random number generator, so that mclapply is reproducible

if (!exists(output)) dir.create(output)
outdir <- normalizePath(output, mustWork = F)

setwd(outdir)

##For now, we are using the regions to exclude and the regions to focus on extracted from the original package.
data(exclude_regions) ##This contains three tricky genes/regions in terms of CNV/CNA. HLA (obvious reasons), CEACAM-PSG (3 to 80% copy number human polymorphism), and GSTT1, with also know CNV polymorphism in humans pops).
data(detail_regions) ##This contains classical cancer genes. I should expand to add BE fragile sites, for example.

##Default minimum 15 probes per bin, minimum bin size 50kb and maximum 5M. Bins are not only of min_size because the probes may be more sparse.
cli_alert_info("Running CNV analysis")
annotations = CNV.create_anno(
  array_type = "EPIC",
  exclude_regions = exclude_regions,
  detail_regions = detail_regions
)

##List of control samples
isAControl =  tolower(pData(preprocessedData)$Group) %in% c("normal", "control")

##List of problem samples
samplesForAnalysis = pData(preprocessedData)$Sample_Name[!isAControl]

##Intensity data
combinedIntensitiesFun = CNV.load(preprocessedDataMS)

##Weird problem with some annotations that do not overlap with the array data. Found the solution here https://support.bioconductor.org/p/97909/
validAnnotations = annotations
validAnnotations@probes = subsetByOverlaps(validAnnotations@probes, granges(preprocessedData))

#Freeing some memory
rm(preprocessedData)
rm(annotations)
#Function to generate the results
getTheFit = function (x, data, annotations, controlMask) {
  sample = x
  out = CNV.fit(data[sample], data[controlMask], annotations) ##Calculates relative log2ratio using the normal
  out = CNV.segment(CNV.detail(CNV.bin(out)))
  out
}

##Main loop to generate results
results = mclapply(
  samplesForAnalysis,
  mc.cores = mc.cores,
  FUN = getTheFit,
  data = combinedIntensitiesFun,
  annotations = validAnnotations,
  controlMask = isAControl
)

names(results) = samplesForAnalysis

printResults = function (x, outdir) {
  pdf(paste0(outdir,"/",x@name,"_genomePlot.pdf"))
  CNV.genomeplot(x)
  dev.off()
  pdf(paste0(outdir,"/",x@name,"_detailPlot.pdf"))
  CNV.detailplot_wrap(x)
  dev.off()
}

rm(combinedIntensitiesFun)


CNV.call = function (x, alpha = 0.001, nperm = 50000, min.width = 5, undo.splits = "sdundo", undo.SD = 2.2, verbose = 0, min_ploidy = 1.25, max_ploidy = 6, ploidy_step = 0.01, min_cellularity = 0.2, max_cellularity = 1, cellularity_step = 0.01, distance_function = "MAD", ...) {
  if (class(x) != "CNV.analysis" || length(x@fit) == 0)
    stop("fit unavailable, run CNV.fit")
  
  # Rascal does not like log2 transformed values, but conumee generates them by default, centered by an optimized shift value
  # We are unshifting and undoing the log2 transformation when building the input for DNAcopy
  
  #Since DNAcopy expects log2 ratios, we are using them to segment the data, but then use the actual relative copy numbers for rascal
  DNAcopyInput <- DNAcopy::CNA(
      genomdat = x@bin$ratio[names(x@anno@bins)],
      chrom = as.vector(seqnames(x@anno@bins)),
      maploc = values(x@anno@bins)$midpoint,
      data.type = "logratio",
      sampleid = x@name
    )
  
  ## Segmentation using DNAcopy exactly as conumee does
  DNAcopyOutput <- DNAcopy::segment(
      x = DNAcopyInput,
      alpha = alpha ,
      nperm = nperm,
      min.width = min.width,
      undo.splits = undo.splits,
      undo.SD = undo.SD,
      verbose = verbose,
      ...
    )
  
  ## DNAcopy::segment returns results by segment, but we want to access the by-bin data to calculate Rascal's weights and to use relative copy numbers isntead of log2 rates
  ## Doing some data shuffling to generate this info
  DNAcopyOutputDT = data.table(DNAcopyOutput$output)
  DNAcopyOutputDT = DNAcopyOutputDT[mixedorder(chrom),]
  
  ##Mixedsort, could still not be enough for some genomes, so we better make sure the order is the same before merging the data
  if (!identical(as.character(unique(decode(
    seqnames(x@anno@bins)
  ))), unique(DNAcopyOutputDT[, chrom]))) {
    stop(
      "Sorting problems between conumee and DNAcopy. This is a bug that needs to be fixed here. Contact the author at dmalload@asu.edu if you are not him."
    )
  }
  
  dataForRascalByBin = data.table(
    chromosome = gsub("chr", "", decode(seqnames(x@anno@bins))),
    start = start(x@anno@bins),
    end = end(x@anno@bins),
    shiftedlogRatio = x@bin$ratio[names(x@anno@bins)] -
      x@bin$shift,
    relativeCopyNumber = 2 ** x@bin$ratio[names(x@anno@bins)],
    segmentedlogRatio = unlist(sapply(
      seq(1, length(DNAcopyOutputDT$ID)),
      FUN = function(x, data) {
        rep(data$seg.mean[x], data$num.mark[x])
      },
      data = DNAcopyOutputDT
    )),
    segmentID = unlist(sapply(
      seq(1, length(DNAcopyOutputDT$ID)),
      FUN = function(x, data) {
        rep(x, data$num.mark[x])
      },
      data = DNAcopyOutputDT
    ))
  )
  
  ##SegmentedLogRatio must be equal to meanLogRatio or otherwise something wrong happened
  
  dataForRascalByBin[, `:=`(length = end - start + 1)]
  
  dataForRascalBySegment = dataForRascalByBin[, .(
    chromosome = first(chromosome),
    start = min(start),
    end = max(end),
    medianRelativeCopyNumber =
      median(relativeCopyNumber),
    meanRelativeCopyNumber =
      mean(relativeCopyNumber),
    meanLogRatio = mean(shiftedlogRatio),
    medianLogRatio = median(shiftedlogRatio),
    length = sum(length),
    weight = sum(length) / median(length),
    nbins = .N
  ),
  by = segmentID]
  
  #Reformatting the previous to match Rascal's format
  #Median normalization (this is important and Rascal does it also in the byBin copy number estimates for plotting (normalizing by the mean segment))
  finalRascalFormat = dataForRascalBySegment[, .(
    chromosome,
    start,
    end,
    copy_number = meanRelativeCopyNumber /
      median(dataForRascalBySegment$meanRelativeCopyNumber),
    bin_count = nbins,
    sum_of_bin_lengths = length,
    weight
  )]
  
  
  solutions <- find_best_fit_solutions(
    finalRascalFormat$copy_number,
    finalRascalFormat$weight,
    min_ploidy = min_ploidy,
    max_ploidy = max_ploidy,
    ploidy_step = ploidy_step,
    min_cellularity = min_cellularity,
    max_cellularity = max_cellularity,
    cellularity_step = cellularity_step,
    distance_function = distance_function
  )
  
  return(
    list(
      dataBySegment = dataForRascalBySegment,
      dataByBin = dataForRascalByBin,
      solutions = solutions
    )
  )
}

CNV.selectSolution = function (x, method = c("highCell", "best")) {
  if (length(x$solutions) == 0)
    stop("Solutions unavailable. Run CNV.call")
  
  if (method == "best") {
    solutionN = as.numeric(
      x$solutions %>% rownames_to_column() %>% arrange(distance, desc(cellularity)) %>% slice(1) %>% select(rowname) %>% mutate(rowname =
                                                                                                                                  rowname)
    )
  } else if (method == "highCell") {
    solutionN = as.numeric(
      x$solutions %>% rownames_to_column() %>% arrange(desc(cellularity), distance) %>% slice(1) %>% select(rowname) %>% mutate(rowname =
                                                                                                                                  rowname)
    )
  } else {
    stop(paste0("The strategy ", method, " is not supported"))
  }
  return(solutionN)
}

CNV.getCalls = function (x, n_solution = NA, method = "best", alterationThreshold = 0.25) {
  if (length(x$solutions) == 0 || length(x$dataBySegment) == 0)
    stop("Solutions or by-segment data unavailable. Make sure to run CNV.call first")
  
  if (is.na(n_solution)) {
    warning(
      paste0(
        "Preferred solution not provided. Obtaining the best solution with method ",
        method
      )
    )
    n_solution = CNV.selectSolution(x, method = method)
  }
  
  ploidy = as.numeric(x$solution[n_solution, "ploidy"])
  cellularity = as.numeric(x$solution[n_solution, "cellularity"])
  
  theCalls = x$dataBySegment[, .(chromosome,
                                 start,
                                 end,
                                 relative_copy_number = meanRelativeCopyNumber,
                                 log2_ratio = meanLogRatio)] %>%
    mutate(absolute_copy_number = relative_to_absolute_copy_number(relative_copy_number, ploidy, cellularity)) %>%
    mutate(across(
      c(log2_ratio, relative_copy_number, absolute_copy_number), \(x) round(x, digits = 3)
    )) %>%
    mutate(
      alteration = case_when(
        absolute_copy_number > ploidy + alterationThreshold ~ "Gain",
        absolute_copy_number < ploidy - alterationThreshold ~ "Loss",
        .default = "Normal"
      )
    )
  
  return(list(
    calls = theCalls,
    ploidy = ploidy,
    cellularity = cellularity
  ))
}

finalCalls = mclapply(results, mc.cores = mc.cores, CNV.call)
names(finalCalls) = names(results)

callsTable = rbindlist(lapply(
  names(finalCalls),
  FUN = function(x) {
    return(data.table(
      CNV.getCalls(finalCalls[[x]], method = "highCell")$calls,
      sample = x
    ))
  }
)) 

if (any(callsTable$alteration != "Normal")) {
  mask <- callsTable %>%
    dplyr::filter(!is.na(relative_copy_number) & alteration != "Normal") %>%
    mutate(chromosome = paste0("chr", chromosome)) %>% 
    makeGRangesFromDataFrame # Create a Granges object from the CNVs to intersect it with the probes
  
  blacklist <-
    names(subsetByOverlaps(validAnnotations@probes, mask, ignore.strand = TRUE)) # Get the probes that fall in CNV regions
  
}

cli_alert_info("Printing the results! Just a few more seconds...")

blacklist %>% 
  as.data.frame %>% 
  rename("id" = ".") %>% 
  write.csv(file = paste0(name, ".blacklist.csv"),
            quote = F,
            row.names = F)


write.table(
  callsTable,
  file = paste0(name, "allCNVCalls.csv"),
  quote = F,
  sep = ",",
  row.names = F
)

#Function to generate the graphical output
printResults = function (x, outdir) {
  pdf(paste0(outdir,"/",x@name,"_genomePlot.pdf"))
  CNV.genomeplot(x)
  CNV.detailplot_wrap(x)
  dev.off()
}

if (genome_plot) {
  dir = "cnv_plots"
  if (!exists(dir)) dir.create(dir)
  #Loop to generate the output
  mclapply(results,mc.cores=mc.cores,printResults,outdir=dir)
}

