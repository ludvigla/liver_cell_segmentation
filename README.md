# cell segmentation in HE images

## Installation requirements

The following R packages are required to run the scripts:

`argparse`
`data.table`
`magick`
`magrittr`
`EBImage`
`dplyr`

A simple walkthrough is available in the doc/cell_segmentation.Rmd rmarkdown where 
you can play around with the parameter settings to achieve optimal results.

## Run script from command line

Once the parameter settings have been determines, the cell segmentation workflow can be run
from the terminal:

`Rscript ./scripts/segment_liver_nuclei.R -h`

Example:

`Rscript ./scripts/segment_liver_nuclei.R  --croparea 12700,17800,12000,7800 --x-scale 0.25 --brush-size-thresholding 21 --brush-size-opening 5 --offset-threshold 0.05 --array-type 1k --method 2 /path/to/HE_image.jpg`



![Alt text](results/cells_within_spots.jpeg?raw=true "Title")
