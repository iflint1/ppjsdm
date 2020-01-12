library(usethis)

raw_qld_trees <- read.csv("data-raw/CSIRO_PermanentPlots_Data/data/CSIRO_PermanentPlots_TreeMeasurementData.csv", stringsAsFactors = FALSE)

# There is a typo here, the data in the csv file is 'eP4'
raw_qld_trees$epNumber[74509] <- "ep4"

# A few spatial values are NA, remove them
# which(is.na(raw_csiro_trees$coordinates_x_metres)) to obtain values below
raw_qld_trees <- raw_qld_trees[-c(17759,
                                      47793,
                                      47794,
                                      49796,
                                      49797,
                                      71310,
                                      71311,
                                      71938,
                                      71939,
                                      71940,
                                      74039,
                                      74040,
                                      74041,
                                      74042,
                                      79724,
                                      95325), ]

usethis::use_data(raw_qld_trees, overwrite = TRUE)
