# Load libraries
library(readr)
library(dplyr)
library(geoChronR)
library(jsonlite)

# Load your merged dataset (update the path as needed)
merged_df <- read_csv("merged_df.csv")

# Check the structure
glimpse(merged_df)

# Define archive types with their corresponding age uncertainties
low_uncertainty_archives <- c("coral", "Wood", "GlacierIce")
high_uncertainty_archives <- c("LakeSediment", "MarineSediment", "Speleothem")

bam_results <- list()

for (this_id in unique_ids) {
  
  one_record <- merged_df %>% filter(ID == this_id)
  
  if (nrow(one_record) == 0) next
  
  archive_type <- one_record$archiveType[1]
  
  # Choose ageVar based on archive type
  ageVar <- if (archive_type %in% low_uncertainty_archives) 0.01 else 0.05
  
  year <- one_record$Year
  value <- one_record$Value
  
  n <- length(year)
  p <- 1
  
  param_array <- array(ageVar, dim = c(n, p, 2))
  
  model_list <- list(
    name = "poisson",
    param = param_array,
    ns = 100,  # or whatever ensemble size you prefer
    resize = 0
  )
  
  bam_output <- simulateBam(
    X = matrix(value, ncol = 1),
    t = year,
    model = model_list
  )
  
  Xp <- bam_output$Xp   # Dimensions: [time, ensemble]
  
  # Apply normalization per ensemble member (each column)
  Xp_normalized <- apply(Xp, 2, function(col) {
    col_mean <- mean(col, na.rm = TRUE)
    col_sd <- sd(col, na.rm = TRUE)
    
    # Check for NA mean
    if (is.na(col_mean)) {
      warning("Column has NA mean. Returning NAs.")
      return(rep(NA, length(col)))
    }
    
    # If sd is zero or NA, center only
    if (col_sd == 0 || is.na(col_sd)) {
      warning("Column has zero or NA standard deviation. Centering only.")
      return(col - col_mean)
    }
    
    # Normal case
    return((col - col_mean) / col_sd)
  })
  
  bam_results[[this_id]] <- list(
    ID = this_id,
    archiveType = archive_type,
    bam_output = list(
      tp = bam_output$tp,                
      Xp = Xp,                           
      Xp_normalized = Xp_normalized      
    )
  )
}

# Export to JSON
json_bam <- toJSON(bam_results, pretty = TRUE)

# Save to file (update path if needed)
write(json_bam, "bam_results_normalized.json")
message("Finished exporting bam_results_normalized.json")