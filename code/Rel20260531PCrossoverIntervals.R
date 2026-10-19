setwd("C:/Users/rafae/OneDrive/Documentos/Religion")

suppressPackageStartupMessages({
  library(readr)
})

#### settings you may want to edit
include_national_rows_in_fit <- TRUE

output_file <- "output/pc_crossover_year_intervals.csv"


#### read the pooled parameter-fit table
TBS <- as.data.frame(
  readr::read_csv(
    "input/Religion 20260324 parameter fit.csv",
    show_col_types = FALSE,
    locale = readr::locale(encoding = "UTF-8")
  )
)

names(TBS)[names(TBS) == "CATOLIC"] <- "C"
names(TBS)[names(TBS) == "NONE"] <- "N"
names(TBS)[names(TBS) == "PROTESTANT"] <- "P"


#### read the observed census tables to get the latest national C/P ratios
mx_obs <- as.data.frame(
  readr::read_csv(
    "input/ReligionMexico 20251205.csv",
    show_col_types = FALSE,
    locale = readr::locale(encoding = "UTF-8")
  )
)
names(mx_obs)[1] <- "STATE"
names(mx_obs)[names(mx_obs) == "AGE G"] <- "AGE_G"

br_obs <- as.data.frame(
  readr::read_csv(
    "input/ReligionBrazil 20260311.csv",
    show_col_types = FALSE,
    locale = readr::locale(encoding = "UTF-8")
  )
)
names(br_obs)[1] <- "STATE"
names(br_obs)[names(br_obs) == "AGE G"] <- "AGE_G"


#### pooled country fit, same structure used for the parameter estimates
fit_switching_panel_closed <- function(TBS,
                                       time_col = "t",
                                       state_col = "STATE",
                                       eps = 1e-8) {
  needed <- c("C", "N", "P", time_col, state_col)
  dat <- TBS[, needed]
  dat <- dat[stats::complete.cases(dat), ]

  names(dat)[names(dat) == time_col] <- "time"
  names(dat)[names(dat) == state_col] <- "STATE"

  obs <- stats::aggregate(cbind(C, N, P) ~ STATE + time, data = dat, FUN = mean)

  s <- obs$C + obs$N + obs$P
  obs <- obs[is.finite(s) & s > 0, ]
  s <- obs$C + obs$N + obs$P
  obs$C <- obs$C / s
  obs$N <- obs$N / s
  obs$P <- obs$P / s

  obs$C <- pmax(obs$C, eps)
  obs$N <- pmax(obs$N, eps)
  obs$P <- pmax(obs$P, eps)

  s <- obs$C + obs$N + obs$P
  obs$C <- obs$C / s
  obs$N <- obs$N / s
  obs$P <- obs$P / s

  obs <- obs[order(obs$STATE, obs$time), ]

  ntab <- table(obs$STATE)
  keep_states <- names(ntab[ntab >= 2])
  obs <- obs[obs$STATE %in% keep_states, ]

  if (nrow(obs) == 0) {
    stop("No usable data after cleaning.")
  }

  predict_one_state <- function(df_state, alpha, beta) {
    df_state <- df_state[order(df_state$time), ]
    t0 <- min(df_state$time)

    init <- df_state[df_state$time == t0, c("C", "N", "P")]
    C0 <- init$C[1]
    N0 <- init$N[1]
    P0 <- init$P[1]

    tt <- df_state$time - t0

    a1 <- alpha * tt
    a2 <- beta * tt
    a3 <- rep(0, length(tt))

    m <- pmax(a1, a2, a3)
    e1 <- exp(a1 - m)
    e2 <- exp(a2 - m)
    e3 <- exp(a3 - m)

    denom <- C0 * e1 + N0 * e2 + P0 * e3

    data.frame(
      STATE = df_state$STATE[1],
      time = df_state$time,
      C = C0 * e1 / denom,
      N = N0 * e2 / denom,
      P = P0 * e3 / denom
    )
  }

  objective <- function(par) {
    alpha <- par[1]
    beta <- par[2]

    states <- unique(obs$STATE)
    sse <- 0

    for (st in states) {
      df_state <- obs[obs$STATE == st, ]
      pred <- predict_one_state(df_state, alpha, beta)

      fitdat <- merge(
        df_state,
        pred,
        by = c("STATE", "time"),
        suffixes = c("_obs", "_mod")
      )

      err <- (fitdat$C_obs - fitdat$C_mod)^2 +
        (fitdat$N_obs - fitdat$N_mod)^2 +
        (fitdat$P_obs - fitdat$P_mod)^2

      if (!all(is.finite(err))) {
        return(1e12)
      }

      sse <- sse + sum(err)
    }

    sse
  }

  start_grid <- rbind(
    c(0, 0),
    c(0.25, 0),
    c(0, 0.25),
    c(-0.25, 0),
    c(0, -0.25),
    c(0.5, 0.5),
    c(-0.5, -0.5),
    c(0.5, -0.5),
    c(-0.5, 0.5)
  )

  best_opt <- NULL
  best_val <- Inf

  for (i in seq_len(nrow(start_grid))) {
    opt_i <- try(
      stats::optim(
        par = start_grid[i, ],
        fn = objective,
        method = "L-BFGS-B",
        lower = c(-10, -10),
        upper = c(10, 10)
      ),
      silent = TRUE
    )

    if (!inherits(opt_i, "try-error") &&
        is.finite(opt_i$value) &&
        opt_i$value < best_val) {
      best_opt <- opt_i
      best_val <- opt_i$value
    }
  }

  if (is.null(best_opt)) {
    stop("Optimization failed.")
  }

  alpha_hat <- unname(best_opt$par[1])
  beta_hat <- unname(best_opt$par[2])

  states <- unique(obs$STATE)
  fitted_all <- do.call(
    rbind,
    lapply(states, function(st) {
      predict_one_state(obs[obs$STATE == st, ], alpha_hat, beta_hat)
    })
  )

  fit_table <- merge(
    obs,
    fitted_all,
    by = c("STATE", "time"),
    suffixes = c("_obs", "_mod")
  )

  resid_vec <- c(
    fit_table$C_obs - fit_table$C_mod,
    fit_table$N_obs - fit_table$N_mod,
    fit_table$P_obs - fit_table$P_mod
  )

  sse <- sum(resid_vec^2)
  df_resid <- length(resid_vec) - 2
  sigma2 <- if (df_resid > 0) sse / df_resid else NA_real_

  num_hessian <- function(fn, par, h = 1e-4) {
    k <- length(par)
    H <- matrix(NA_real_, k, k)

    for (i in 1:k) {
      for (j in i:k) {
        ei <- rep(0, k)
        ej <- rep(0, k)
        ei[i] <- 1
        ej[j] <- 1

        fpp <- fn(par + h * ei + h * ej)
        fpm <- fn(par + h * ei - h * ej)
        fmp <- fn(par - h * ei + h * ej)
        fmm <- fn(par - h * ei - h * ej)

        H[i, j] <- (fpp - fpm - fmp + fmm) / (4 * h^2)
        H[j, i] <- H[i, j]
      }
    }

    H
  }

  H <- try(num_hessian(objective, c(alpha_hat, beta_hat)), silent = TRUE)
  vcov_mat <- matrix(NA_real_, 2, 2, dimnames = list(c("alpha", "beta"), c("alpha", "beta")))

  if (!inherits(H, "try-error") && all(is.finite(H)) && is.finite(sigma2)) {
    Hinv <- try(solve(H), silent = TRUE)
    if (!inherits(Hinv, "try-error") && all(is.finite(Hinv))) {
      vcov_mat <- 2 * sigma2 * Hinv
    }
  }

  list(
    coefficients = c(alpha = alpha_hat, beta = beta_hat),
    vcov = vcov_mat
  )
}


#### prepare estimation data
fit_data <- TBS[, c("COUNTRY", "STATE", "t", "C", "N", "P")]
fit_data <- fit_data[stats::complete.cases(fit_data), ]

if (!include_national_rows_in_fit) {
  fit_data <- fit_data[
    !(fit_data$COUNTRY == "MEXICO" & fit_data$STATE == "Estados Unidos Mexicanos") &
      !(fit_data$COUNTRY == "BRAZIL" & fit_data$STATE == "Brasil"),
  ]
}


#### latest observed national ratios
mx_obs <- mx_obs[mx_obs$AGE_G == "Total", ]
mx_base_year <- max(mx_obs$YEAR, na.rm = TRUE)
mx_nat <- mx_obs[mx_obs$STATE == "Estados Unidos Mexicanos" & mx_obs$YEAR == mx_base_year, ]
mx_c <- mx_nat$TOT[mx_nat$RELIGION == "CATOLIC"]
mx_p <- mx_nat$TOT[mx_nat$RELIGION == "PROTESTANT"]
mx_n <- mx_nat$TOT[mx_nat$RELIGION == "NONE"]
mx_total <- mx_c + mx_p + mx_n
mx_c_share <- mx_c / mx_total
mx_p_share <- mx_p / mx_total
mx_ratio <- mx_c / mx_p

br_base_year <- max(br_obs$YEAR, na.rm = TRUE)
br_nat <- br_obs[br_obs$STATE == "Brasil" & br_obs$YEAR == br_base_year, ]
br_c <- br_nat$TOT[br_nat$RELIGION == "CATOLIC"]
br_p <- br_nat$TOT[br_nat$RELIGION == "PROTESTANT"]
br_n <- br_nat$TOT[br_nat$RELIGION == "NONE"]
br_total <- br_c + br_p + br_n
br_c_share <- br_c / br_total
br_p_share <- br_p / br_total
br_ratio <- br_c / br_p


#### fit and compute crossover years
country_specs <- list(
  list(country = "MEXICO", base_year = mx_base_year, cp_ratio = mx_ratio, c_share = mx_c_share, p_share = mx_p_share),
  list(country = "BRAZIL", base_year = br_base_year, cp_ratio = br_ratio, c_share = br_c_share, p_share = br_p_share)
)

results <- list()

for (i in seq_along(country_specs)) {
  spec <- country_specs[[i]]
  message("Fitting ", spec$country, "...")

  df_country <- fit_data[fit_data$COUNTRY == spec$country, ]
  fit <- fit_switching_panel_closed(df_country)

  alpha_hat <- unname(fit$coefficients["alpha"])
  alpha_se <- sqrt(fit$vcov[1, 1])
  alpha_ci_low <- alpha_hat - 1.96 * alpha_se
  alpha_ci_high <- alpha_hat + 1.96 * alpha_se

  cp_ratio <- unname(spec$cp_ratio[1])
  base_year <- unname(spec$base_year[1])
  c_share <- unname(spec$c_share[1])
  p_share <- unname(spec$p_share[1])

  t_cross_est <- -log(cp_ratio) / alpha_hat
  t_cross_low <- -log(cp_ratio) / alpha_ci_low
  t_cross_high <- -log(cp_ratio) / alpha_ci_high

  year_est <- base_year + 10 * t_cross_est
  year_low <- base_year + 10 * t_cross_low
  year_high <- base_year + 10 * t_cross_high

  results[[i]] <- data.frame(
    country = spec$country,
    base_year = base_year,
    catholic_share_base = c_share,
    protestant_share_base = p_share,
    cp_ratio_base = cp_ratio,
    alpha_estimate = alpha_hat,
    alpha_se = alpha_se,
    alpha_ci_low_95 = alpha_ci_low,
    alpha_ci_high_95 = alpha_ci_high,
    crossover_time_decades_estimate = t_cross_est,
    crossover_time_decades_low_95 = t_cross_low,
    crossover_time_decades_high_95 = t_cross_high,
    crossover_year_estimate = year_est,
    crossover_year_low_95 = year_low,
    crossover_year_high_95 = year_high,
    crossover_year_estimate_round = round(year_est),
    crossover_year_low_95_round = round(year_low),
    crossover_year_high_95_round = round(year_high),
    method = "latest observed C/P ratio + alpha hessian interval"
  )
}

results_table <- do.call(rbind, results)
row.names(results_table) <- NULL

write.csv(results_table, output_file, row.names = FALSE, fileEncoding = "Windows-1252")

print(results_table)
message("Wrote: ", output_file)
