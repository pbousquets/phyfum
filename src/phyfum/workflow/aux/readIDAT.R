if (!require("pacman")) install.packages("pacman")
p_load(optparse, cli)

# Define the command line options
option_list <- list(
  make_option(
    c("--input", "-i"),
    type = "character",
    default = NULL,
    help = "Directory containing CSV files (required)"
  ),
  make_option(
    c("--mc.cores", "-c"),
    type = "integer",
    default = 1,
    help = "Number of CPU cores to use [default: %default])"
  ),
  make_option(
    c("--patientInfo", "-p"),
    type = "character",
    default = NULL,
    help = "Patient info file (required,  [default: %default])"
  ),
  make_option(
    c("--output", "-o"),
    type = "character",
    default = "idat_processed",
    help = "Output directory (required,  [default: %default])"
  ),
  make_option(
    c("--name", "-n"),
    type = "character",
    default = "myexperiment",
    help = "Experiment name (required,  [default: %default])"
  )
)

# Parse the command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Libraries
cli_alert_info("Preparing the environment")

p_load(minfi)

# Check if libraryd arguments are provided
if (is.null(opt$input)) {
  cat("Error: Please provide an input directory.\n")
  cat("Run with --help for usage information.\n")
  quit(status = 1)
}

input_dir <- opt$input
mc.cores <- opt$mc.cores
patientInfo <- opt$patientInfo
output <- opt$output
name <- opt$name
rm(opt)

if (!exists(output)) dir.create(output)
outdir <- normalizePath(output, mustWork = F)

cli_alert_info("Reading the data")

# Read raw data
setwd(input_dir)

targets = read.csv(patientInfo, header = T)
rawData = read.metharray.exp(targets = targets, force = T, recursive = T) #Activate force in case there are samples with of different arrays

setwd(outdir)
cli_alert_info("Normalizing")

##Pre-processing data. There are multiple options here, Illumina, SWAN, Quantile, Noob, Funnorm. We use Noob for now. Funnorm seems to be clearly the best
preprocessedData = preprocessFunnorm(rawData, ratioConvert = FALSE)

 ##preprocessedData needs to be converted from a genomicMethylSet back to a MethylSet to use CNV.load.
preprocessedDataMS = MethylSet(
  Meth = getMeth(preprocessedData),
  Unmeth = getUnmeth(preprocessedData),
  colData = colData(preprocessedData),
  preprocessMethod = preprocessMethod(preprocessedData),
  annotation = annotation(preprocessedData),
  metadata = metadata(preprocessedData)
)

cli_alert_info("Writing the final output!")

write.table(getBeta(preprocessedDataMS), file = paste0(name, ".betas.csv"),
            quote = F,
            sep = ",",
            row.names = T)

write.table(getMeth(preprocessedDataMS), file = paste0(name, ".M.csv"),
            quote = F,
            sep = ",",
            row.names = T)

write.table(getUnmeth(preprocessedDataMS), file = paste0(name, ".U.csv"),
            quote = F,
            sep = ",",
            row.names = T)

save.image(paste0(name, ".RData"),)