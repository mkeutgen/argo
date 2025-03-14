# Install if not already installed:
# install.packages("fitdistrplus")
# install.packages("poweRlaw")
# install.packages("evir")
# install.packages("tidyverse")

library(fitdistrplus)
library(poweRlaw)
library(evir)       # for GPD
library(tidyverse)


get_fitted_pdf_data <- function(sub_df, n_points = 500, tail_prob = 0.8) {
  # sub_df is the subset for one region
  
  # 1) Extract the numeric vector
  x_vals <- sub_df$PRES_ADJUSTED
  
  # 2) Fit standard distributions via fitdistrplus
  fit_exp       <- fitdist(x_vals, "exp")
  fit_gamma     <- fitdist(x_vals, "gamma")
  fit_lognormal <- fitdist(x_vals, "lnorm")
  fit_cauchy    <- fitdist(x_vals, "cauchy",
                           start = list(location = mean(x_vals), 
                                        scale    = sd(x_vals)))
  
  # 3) Fit a power law using poweRlaw
  #    Decide if your data is "continuous" or "discrete" for poweRlaw:
  #    - If your data are large integer values or truly continuous, use conpl$new().
  #    - If your data are counts/small discrete, use displ$new().
  pl_model <- conpl$new(x_vals)       # or displ$new(x_vals), if appropriate
  est_xmin <- estimate_xmin(pl_model)
  pl_model$setXmin(est_xmin)
  est_pl   <- estimate_pars(pl_model)
  pl_model$setPars(est_pl)
  
  alpha_hat <- pl_model$getPars()
  xmin_hat  <- pl_model$getXmin()
  
  # For a continuous power law, PDF for x >= xmin is:
  powerlaw_pdf <- function(x, alpha, xmin) {
    ifelse(x < xmin, 0, alpha * xmin^alpha * x^(-alpha - 1))
  }
  

  # 5) Create a sequence of x-values spanning the data range
  x_seq <- seq(
    from = max(0, min(x_vals)),  # assume 0 as lower bound if needed
    to   = max(x_vals),
    length.out = n_points
  )
  
  # 6) Compute PDFs from each fitted model
  exp_density       <- dexp(x_seq, rate = fit_exp$estimate["rate"])
  gamma_density     <- dgamma(x_seq,
                              shape = fit_gamma$estimate["shape"],
                              rate  = fit_gamma$estimate["rate"])
  lognormal_density <- dlnorm(x_seq,
                              meanlog = fit_lognormal$estimate["meanlog"],
                              sdlog   = fit_lognormal$estimate["sdlog"])
  cauchy_density    <- dcauchy(x_seq,
                               location = fit_cauchy$estimate["location"],
                               scale    = fit_cauchy$estimate["scale"])
  
  # For the power law: shift by + xmin if needed, or define directly:
  pl_density <- powerlaw_pdf(x_seq, alpha_hat, xmin_hat)
  

  # 7) Create a long-format data frame of x vs. PDF for each distribution
  density_df <- tibble(
    x         = x_seq,
    Exp       = exp_density,
    Gamma     = gamma_density,
    LogNormal = lognormal_density,
    Cauchy    = cauchy_density,
    PL        = pl_density,
  ) %>%
    pivot_longer(
      cols      = c("Exp", "Gamma", "LogNormal", "Cauchy", "PL"),
      names_to  = "Distribution",
      values_to = "Density"
    )
  
  # 8) Compute GOF measures using gofstat() for the fitdist-based models
  gof_exp       <- gofstat(fit_exp)
  gof_gamma     <- gofstat(fit_gamma)
  gof_lognormal <- gofstat(fit_lognormal)
  gof_cauchy    <- gofstat(fit_cauchy)
  
  # For the power law, we can’t use gofstat() directly. We set NA placeholders.
  # You could implement custom tests or the poweRlaw bootstrap KS if desired.
  
  gof_df <- data.frame(
    Distribution = c("Exp", "Gamma", "LogNormal", "Cauchy", "PL"),
    AIC          = c(gof_exp$aic, gof_gamma$aic, gof_lognormal$aic, 
                     gof_cauchy$aic, NA),
    BIC          = c(gof_exp$bic, gof_gamma$bic, gof_lognormal$bic, 
                     gof_cauchy$bic, NA),
    KS           = c(gof_exp$ks,  gof_gamma$ks,  gof_lognormal$ks,  
                     gof_cauchy$ks,  NA),
    CvM          = c(gof_exp$cvm, gof_gamma$cvm, gof_lognormal$cvm, 
                     gof_cauchy$cvm, NA),
    AD           = c(gof_exp$ad,  gof_gamma$ad,  gof_lognormal$ad,  
                     gof_cauchy$ad, NA)
  )
  
  # Optionally store some key parameters in the density data:
  density_df$pl_xmin  <- xmin_hat
  density_df$pl_alpha <- alpha_hat

  # Return both the density data and GOF measures
  return(list(density_df = density_df, gof_df = gof_df))
}

fitted_list <- pres_data %>%
  group_by(region) %>%
  do( fit_out = get_fitted_pdf_data(.) )

# "fit_out" is now a list column with two sub-elements: density_df and gof_df.
# We can "unnest" them:

fitted_list2 <- fitted_list %>% 
  tidyr::unnest_wider(fit_out)  # splits into "density_df" and "gof_df" columns

# Extract the combined density data for plotting:
fitted_df <- fitted_list2 %>% 
  unnest(cols = density_df)

# Extract the GOF summary for each region:
gof_table <- fitted_list2 %>%
  unnest(cols = gof_df)

# Optional: Inspect the GOF table by region
gof_table %>% 
  arrange(region, Distribution)


ggplot() +
  geom_histogram(
    data = pres_data,
    aes(x = PRES_ADJUSTED, y = ..density.., fill = region),
    bins = 30, color = "black", fill = "white", alpha = 0.8
  ) +
  geom_line(
    data = fitted_df,
    aes(x = x, y = Density, color = Distribution),
    size = 1
  ) +
  facet_wrap(~ region, scales = "free_y") +
  scale_color_manual(values = c(
    "Exp"       = "red",
    "Gamma"     = "blue",
    "LogNormal" = "orange",
    "Cauchy"    = "purple",
    "PL"        = "darkgreen",  # Add colors for new distributions
    "GPD"       = "brown"
  )) +
  labs(title = "Fitted Distributions by Region",
       x     = "Pressure Adjusted",
       y     = "Density") +
  theme_minimal()


# Plot each region's data + fitted CDF
cdf_plot <- ggplot() +
  # Empirical CDF, one for each region
  stat_ecdf(
    data = pres_data,
    aes(x = PRES_ADJUSTED,color=region),
    size = 1,
    alpha = 0.5
  )+geom_hline(yintercept = 0.80)
