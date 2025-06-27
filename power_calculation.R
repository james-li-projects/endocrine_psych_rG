########################################
# Simplified power calculation function with fixed parameters
calc_power_fixed <- function(n1, n2) {
  # Fixed parameters
  hsq1 <- 0.1
  hsq2 <- 0.1
  rg <- 0.15
  rp <- 0.15
  overlap <- FALSE
  alpha <- 0.05
  var_pi <- 2e-5
  
  # Compute variance of rg based on overlap status
  var_rg <- if (overlap) {
    ((1 - rg * rp)^2 + (rg - rp)^2) / (hsq1 * hsq2 * n1^2 * var_pi)
  } else {
    (rg^2 * (n1^2 * hsq1^2 + n2^2 * hsq2^2) + 2 * hsq1 * hsq2 * n1 * n2) /
      (2 * hsq1^2 * hsq2^2 * n1^2 * n2^2 * var_pi)
  }
  
  # Non-centrality parameter
  ncp <- rg^2 / var_rg
  
  # Compute power
  power <- pchisq(qchisq(alpha, df = 1, lower.tail = FALSE), df = 1, ncp = ncp, lower.tail = FALSE)
  
  return(power)
}


#######################################
# Psychiatric traits (x-axis labels and sample sizes)
x_labels <- c(
  "Anorexia Nervosa",
  "Tourette's Syndrome",
  "Obsessive-Compulsive symptoms",
  "Schizophrenia",
  "Attention-Deficit Hyperactivity Disorder",
  "Bipolar Disorder (all)",
  "Type 1 Bipolar Disorder",
  "Type 2 Bipolar Disorder",
  "Major Depressive Disorder",
  "Panic Disorder"
)

x_vals <- c(72517, 14307, 33943, 320404, 225534, 413466, 396609, 378330, 807553, 10240)

# Hormone traits (y-axis labels and sample sizes)
y_labels <- c(
  "Sex Hormone-Binding Globulin",
  "Estradiol",
  "Testosterone",
  "TSH",
  "FT4",
  "FT3",
  "TT3",
  "FT3/FT4",
  "TT3/FT4"
)

y_vals <- c(381526, 67623, 381081, 271040, 119120, 59061, 15829, 51095, 15510)

# Define power calculation function (unchanged)
calc_power_fixed <- function(n1, n2) {
  hsq1 <- 0.1
  hsq2 <- 0.1
  rg <- 0.1
  rp <- 0.1
  overlap <- FALSE
  alpha <- 0.05
  var_pi <- 2e-5
  
  if (overlap) {
    var_rg <- ((1 - rg * rp)^2 + (rg - rp)^2) / (hsq1 * hsq2 * n1^2 * var_pi)
  } else {
    var_rg <- (rg^2 * (n1^2 * hsq1^2 + n2^2 * hsq2^2) + 2 * hsq1 * hsq2 * n1 * n2) / 
      (2 * hsq1^2 * hsq2^2 * n1^2 * n2^2 * var_pi)
  }
  ncp <- rg^2 / var_rg
  power <- pchisq(qchisq(alpha, df = 1, lower.tail = FALSE), df = 1, ncp = ncp, lower.tail = FALSE)
  return(power)
}

# Create all combinations
power_DF <- expand.grid(x_index = seq_along(x_labels), y_index = seq_along(y_labels))

# Add sample sizes and labels
power_DF$x_label <- x_labels[power_DF$x_index]
power_DF$x_n <- x_vals[power_DF$x_index]
power_DF$y_label <- y_labels[power_DF$y_index]
power_DF$y_n <- y_vals[power_DF$y_index]

# Calculate power for each pair
power_DF$power <- mapply(calc_power_fixed, power_DF$x_n, power_DF$y_n)

# Drop index columns
power_DF <- power_DF[, c("x_label", "x_n", "y_label", "y_n", "power")]

# View a preview
head(power_DF)





#############################################
# Hormone labels and sample sizes
hormone_labels <- c(
  "Sex Hormone-Binding Globulin",
  "Estradiol",
  "Testosterone",
  "TSH",
  "FT4",
  "FT3",
  "TT3",
  "FT3/FT4",
  "TT3/FT4"
)

hormone_n <- c(381526, 67623, 381081, 271040, 119120, 59061, 15829, 51095, 15510)

# Create all pairwise combinations (including diagonal)
hormone_pairs <- expand.grid(
  h1_index = seq_along(hormone_labels),
  h2_index = seq_along(hormone_labels)
)

# Map sample sizes and labels
hormone_pairs$h1_label <- hormone_labels[hormone_pairs$h1_index]
hormone_pairs$h1_n <- hormone_n[hormone_pairs$h1_index]
hormone_pairs$h2_label <- hormone_labels[hormone_pairs$h2_index]
hormone_pairs$h2_n <- hormone_n[hormone_pairs$h2_index]

# Calculate power for each hormone pair
hormone_pairs$power <- mapply(calc_power_fixed, hormone_pairs$h1_n, hormone_pairs$h2_n)

# Keep relevant columns
power_DF_only_hormone <- hormone_pairs[, c("h1_label", "h1_n", "h2_label", "h2_n", "power")]

# Preview
head(power_DF_only_hormone)





#################################

library(dplyr)

power_DF_only_hormone_unique <- power_DF_only_hormone %>%
  filter(h1_label != h2_label) %>%
  rowwise() %>%
  mutate(
    h1 = min(h1_label, h2_label),
    h2 = max(h1_label, h2_label),
    n1 = ifelse(h1_label < h2_label, h1_n, h2_n),
    n2 = ifelse(h1_label < h2_label, h2_n, h1_n)
  ) %>%
  ungroup() %>%
  distinct(h1, n1, h2, n2, .keep_all = TRUE) %>%
  select(h1_label = h1, h1_n = n1, h2_label = h2, h2_n = n2, power)

sum(power_DF_only_hormone_unique$power>=0.8)
nrow(power_DF_only_hormone_unique)


sum(power_DF$power>=0.8)
nrow(power_DF)


