#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library("argparse")
  library("data.table")
  library("magick")
  library("magrittr")
  library("EBImage")
  library("ggplot2")
  library("dplyr")
})

# create parser object
docstring <- "Takes an H&E image as input and runs cell segmentation"
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-v", "--verbose", action = "store_true", default = TRUE,
                    help = "Print extra output [default]")
parser$add_argument("-o", "--outdir", default = NULL,
                    help = "Output directory name")
parser$add_argument("-s", "--spotfile", default = NULL,
                    help = "Path to spotfile")
parser$add_argument("--x-scale", default = 0.25, type = "double",
                    help = "Scaling factor image width [numeric]. Need to be between 0 and 1.")
parser$add_argument("--croparea", default = NULL,
                    help = "Rectangular crop area defined by width*height+x_offset+y_offset [string]. Example format: 1000,1000,200,200")
parser$add_argument("--method", default = "1", type = "character",
                    help = "Segmentation method (choices: '1', '2')")
parser$add_argument("--offset-threshold", default = 0.05, type = "double",
                    help = "Offset for thresholding [numeric]")
parser$add_argument("--brush-size-thresholding", default = 5, type = "integer",
                    help = "Brush size for thresholding [numeric]")
parser$add_argument("--brush-size-opening", default = 3, type = "integer",
                    help = "Brush size for opening [numeric]")
parser$add_argument("--cell-threshold-min", default = 15, type = "integer",
                    help = "Minimum cell size [numeric]")
parser$add_argument("--cell-threshold-max", default = 500, type = "integer",
                    help = "Maximum cell size [numeric]")
parser$add_argument("--array-type", default = "2k", type = "character",
                    help = "Set array type ['1k' or '2k']")
parser$add_argument("--skip-seg", action = "store_true", default = FALSE,
                    help = "Skip segmentation process")
parser$add_argument("--disable-watershed", action = "store_false", default = TRUE,
                    help = "Should watershedding be used? This will slow down the segmentation significantly, but increases the performance [default]")
parser$add_argument("--spotfile-scale-factor", default = 0.1, type = "double",
                    help = "Scale factor used to map coordinates from HE image to spotfile coordinate system [numeric]")
parser$add_argument("--export-intermediate-files", action = "store_true", default = FALSE,
                    help = "Export image of segmented nuclei and nuclei outlined on H&E image [default]")

parser$add_argument("image", nargs = 1, help = "H&E image")

# Add arguments
args <- parser$parse_args()
image <- args$image

# Set verbosity
verbose <- args$verbose

# Set output directory
outdir <- ifelse(is.null(args$outdir), "./", args$outdir)
if (!dir.exists(outdir)) {
  dir.create(outdir)
}

# Define spot radius
fct <- ifelse(args$array_type == "1k", 0.25, 50/(sqrt(2)*100))

if (verbose) {
  cat("\n")
  cat("image =", image, "\n")
  cat("outdir =", outdir, "\n")
  cat("spotfile =", args$spotfile, "\n")
  cat("spotfile scaling factor =", args$spotfile_scale_factor, "\n")
  cat("x scaling factor =", args$x_scale, "\n")
  cat("croparea =", args$croparea, "\n")
  cat("segmentation method =", args$method, "\n")
  cat("offset for thresholding =", args$offset_threshold, "\n")
  cat("brush size thresholding =", args$brush_size_thresholding, "\n")
  cat("brush size opening =", args$brush_size_opening, "\n")
  cat("cell threshold min =", args$cell_threshold_min, "\n")
  cat("cell threshold max =", args$cell_threshold_max, "\n")
  cat("array type =", args$array_type, "\n")
  cat("spot radius =", fct, "spot-spot distances\n")
  cat("watershed =", args$disable_watershed, "\n")
  cat("export intermediate files =", args$export_intermediate_files, "\n\n")
}

if (verbose) {
  cat("Reading", image, "...\n")
}
im <- image_read(path = image)
iminf <- image_info(im)
# read crop rectangle
croprect <- args$croparea
# Check scaling factor
stopifnot(args$x_scale < 1 & args$x_scale > 0)
if (!is.null(croprect)) {
  if (verbose) {
    cat("Cropping image to a rectangle specified by:", croprect, "...\n")
  }
  rect <- as.numeric(unlist(strsplit(croprect, ",")))
  stopifnot(length(rect) == 4)
  w <- rect[1]; h <- rect[2]; xoff <- rect[3]; yoff <- rect[4]
  im <- im %>%
    image_crop(geometry = geometry_area(width = w, height = h, x_off = xoff, y_off = yoff))
  iminf <- image_info(im)
  print(paste0("Image width: ", iminf$width))
  print(paste0("Image width after scaling: ", args$x_scale*iminf$width))
  im <- im %>%
    image_scale(paste0(round(args$x_scale*iminf$width))) %>%
    as_EBImage()
} else {
  im <- im %>%
    image_scale(paste0(args$x_scale)) %>%
    as_EBImage()
}

if (!args$skip_seg) {
  
  # Clean cells function
  clean_cells <- function (
    imthreshold,
    thr,
    verbose = FALSE
  ) {
    if (class(imthreshold) != "Image") stop(paste0("Invalid input format: ", class(imthreshold)))
    if (verbose) cat("Cleaning up unwanted speckles ... \n")
    imthreshold <- bwlabel(imthreshold)
    areas <- table(imthreshold)[-1]
    inds <- which(areas < thr[1] | areas > thr[2])
    imclean <- rmObjects(x = imthreshold, index = inds)
    return(imclean)
  }
  
  # Watershed function
  watershed_cells <- function (
    imclean,
    tol = 0.1,
    verbose = FALSE
  ) {
    if (class(imclean) != "Image") stop(paste0("Invalid input format: ", class(imclean)))
    if (verbose) cat("Applying watershed ... \n")
    imwatershed <- EBImage::watershed(x = EBImage::distmap(imclean), tolerance = tol)
    return(imwatershed)
  }
  
  
  m1 <- function(imgray, args) {
    nmask = thresh(imgray, w = args$brush_size_thresholding, h = args$brush_size_thresholding, offset = args$offset_threshold)
    nmask = opening(nmask, makeBrush(args$brush_size_opening, shape = 'disc'))
    nmask = fillHull(nmask)
    return(nmask)
  }
  
  m2 <- function(imgray, args) {
    f = makeBrush(args$brush_size_thresholding, shape='disc', step=FALSE)
    f = f/sum(f)
    nmask <- imgray > (filter2(imgray, f, boundary="replicate") + args$offset_threshold)
    nmask = opening(nmask, makeBrush(args$brush_size_opening, shape='disc'))
    nmask = fillHull(nmask)
    return(nmask)
  }
  
  # Define segmentation function
  seg_cells <- function(
    im, 
    args
  ) {
    
    imrb <- rgbImage(red = channel(im, "red"), blue = channel(im, "blue"))
    imgray <- channel(imrb, "gray")
    imgray <- 1 - imgray
    imgray <- imgray^2
    
    if (args$method == "1") {
      nmask <- m1(imgray, args)
    } else if (args$method == "2") {
      nmask <- m2(imgray, args)
    }
    
    nmask = fillHull(nmask)
    cells_segmented <- bwlabel(nmask)
    cells_cleaned <- clean_cells(imthreshold = cells_segmented, thr = c(args$cell_threshold_min, args$cell_threshold_max), verbose = TRUE)
    if (args$disable_watershed) {
      cells_watershed <- watershed_cells(imclean = cells_cleaned, verbose = TRUE)
    } else {
      cells_watershed <- cells_cleaned
    }
    return(cells_watershed)
  }
  
  # Run workflow
  if (verbose) {
    cat("Running segmentation workflow ...\n")
  }
  cleaned_cells_section <- seg_cells(im, args)
  
  # Export segmented cells object
  if (verbose) {
    cat("Export segmented cells ...\n")
  }
  saveRDS(object = cleaned_cells_section, file = paste0(outdir, "/segmented_cells.Rds"))
  
  
  # Export images
  if (args$export_intermediate_files) {
    cells_colored_section <- paintObjects(cleaned_cells_section, im, col = "#FFA500")
    writeImage(x = cells_colored_section, type = "jpeg", files = paste0(outdir, "/outlined_cells.jpeg"), quality = 100)
    writeImage(x = colorLabels(cleaned_cells_section), type = "jpeg", files = paste0(outdir, "/labelled_cells.jpeg"), quality = 100)
  }
  if (verbose) {
    cat("Segmentation complete!\n\n")
  }
} else {
  cat("\n[Skipped segmentation!!!] Reading segmented output to extract numbers of cells per spot...\n\n")
  cleaned_cells_section <- readRDS(file = paste0(outdir, "/segmented_cells.Rds"))
}


if (!is.null(args$spotfile)) {
  if (verbose) {
    cat("Computing cell features ...\n")
  }
  # Extract centroid coordinates for cells
  fts.moment <- computeFeatures.moment(cleaned_cells_section)
  
  # Create a cell coordinate table 
  cell_coordinates <- data.frame(fts.moment[, c("m.cx", "m.cy")])
  
  original_HE_width <- w # width of HE used for ST spot detector (pixels)
  downscaled_HE_width <- round(args$x_scale*iminf$width) # width of downscaled HE image used for this segmentation process
  sf <- original_HE_width/downscaled_HE_width # Scalefactor
  
  # Convert cell coordinates to fit to the same coordinate system as the spots
  cell_coordinates$m.cx <- ((cell_coordinates$m.cx*sf) + xoff)*args$spotfile_scale_factor
  cell_coordinates$m.cy <- ((cell_coordinates$m.cy*sf) + yoff)*args$spotfile_scale_factor
  
  # Spot coordinates
  if (verbose) {
    cat("Loading spot coordinates...\n")
  }
  spot_coordinates <- read.table(file = args$spotfile, header = T, sep = "\t")
  
  # Minimum distances
  if (verbose) {
    cat("Calculating average spot center-to-center distance ...\n")
  }
  distMat <- dist(spot_coordinates[, c("pixel_x", "pixel_y")])
  distMat <- as.matrix(distMat)
  diag(distMat) <- Inf
  center_to_center_distance <- mean(apply(distMat, 2, min))
  
  if (verbose) {
    cat("Defining spot radius...\n")
  }
  spot_radius <- center_to_center_distance*fct
  
  # Export H&E with drawn spots
  if (args$export_intermediate_files) {
    
    # Define circle draw function
    plotCircle <- function(x, y, r) {
      angles <- seq(0, 2*pi, length.out = 360)
      lines(r*cos(angles) + x, r*sin(angles) + y)
    }
    
    # draw circles on H&E image and export
    jpeg(filename = paste0(outdir, "/HE_spots.jpeg"), width = dim(im)[1], height = dim(im)[2], res = 300)
    display(im, method = "raster")
    apply(spot_coordinates[, c("pixel_x", "pixel_y")], 1, function(df) {
      plotCircle(x = ((df[1] - xoff*args$spotfile_scale_factor)/sf)/args$spotfile_scale_factor, 
                 y = ((df[2] - yoff*args$spotfile_scale_factor)/sf)/args$spotfile_scale_factor, 
                 r = (spot_radius/sf)/args$spotfile_scale_factor)
    })
    dev.off()
  }
  
  # Calculate pairwise distances
  if (verbose) {
    cat("Calculating spot-to-cell distances ...\n")
  }
  set1 <- spot_coordinates[, c("pixel_x", "pixel_y")]
  set2 <- cell_coordinates
  mindist <- apply(set1, 1, function(x) {
    sqrt(colSums((t(set2[, 1:2]) - x)^2))
  })
  spotids <- paste0(spot_coordinates$x, "x", spot_coordinates$y)
  
  arrind <- which(mindist <= spot_radius, arr.ind = TRUE)
  df <- data.frame(dist = apply(arrind, 1, function(x) {mindist[x[1], x[2]]}),
                   ID = arrind[, "row"], spot = spotids[arrind[, "col"]])
  cell_coordinates$within <- "no"
  cell_coordinates$within[df$ID] <- "yes"
  cell_coordinates$spot <- NA
  cell_coordinates$spot[df$ID] <- df$spot
  
  if (args$export_intermediate_files) {
    p <- ggplot(cell_coordinates, aes(m.cx, m.cy, color = within)) +
      geom_point(size = 1) +
      theme_void() +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_reverse(expand = c(0, 0))
    
    jpeg(filename = paste0(outdir, "/cells_within_spots.jpeg"), width = dim(im)[1], height = dim(im)[2], res = 300)
    print(p)
    dev.off()
    
    all_inds <- 1:length(table(cleaned_cells_section))
    rm_inds <- setdiff(all_inds, df$ID)
    cleaned_cells_section_filtered <- rmObjects(cleaned_cells_section, index = rm_inds, reenumerate = FALSE)
    writeImage(x = cleaned_cells_section_filtered, type = "jpeg", files = paste0(outdir, "/cells_on_spots.jpeg"), quality = 100)
  }
  
  # Cells per group
  if (verbose) {
    cat("Summarizing cells per spot ...\n")
  }
  cell_coordinates_subset <- na.omit(cell_coordinates)
  cells_per_spot <- cell_coordinates_subset %>% 
    dplyr::group_by(within, spot) %>%
    dplyr::summarize(count = n())
  
  write.table(x = cells_per_spot, file = paste0(outdir, "/cells_per_spot_table.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)
  if (verbose) {
    cat("Pipeline finished!\n")
  }
}