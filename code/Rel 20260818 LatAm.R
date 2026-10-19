setwd("C:/Users/rafae/OneDrive/Documentos/Religion")

suppressPackageStartupMessages({
  library(readr)
})

### This script reads the Latin America template, uses the observed 2014 religion
### shares as the baseline for each country, and fills the reserved forecast rows
### from 2015 to 2050.

### The lower bound applies the Mexico switching parameters.
### The upper bound applies the Brazil switching parameters.
### For Mexico and Brazil themselves, both bounds use their own country estimates.
### The mean is the simple average between the lower and upper bounds.
### The LatAm rows are computed as population-weighted aggregates of the country rows.

input_file <- "input/ReligionLatAm20260818.csv"
output_file <- "outputs/ReligionLatAmForecasts.csv"
mirror_output_file <- "output/ReligionLatAmForecasts.csv"

dir.create("outputs", showWarnings = FALSE)

db <- as.data.frame(
  readr::read_csv(
    input_file,
    show_col_types = FALSE,
    locale = readr::locale(encoding = "UTF-8")
  )
)

### Parameters with Unaffiliated as the reference group.
### These come from the previously estimated country models.
rho_mexico_c <- -0.183255220747426
rho_mexico_p <-  0.067765122086337

rho_brazil_c <- -0.495485629074094
rho_brazil_p <-  0.0589895036188038


### Exact national census anchors for Mexico and Brazil.
mexico_obs <- as.data.frame(
  readr::read_csv(
    "input/ReligionMexico 20251205.csv",
    show_col_types = FALSE,
    locale = readr::locale(encoding = "UTF-8")
  )
)
brazil_obs <- as.data.frame(
  readr::read_csv(
    "input/ReligionBrazil 20260311.csv",
    show_col_types = FALSE,
    locale = readr::locale(encoding = "UTF-8")
  )
)

names(mexico_obs)[1] <- "STATE"
names(brazil_obs)[1] <- "STATE"
names(mexico_obs)[names(mexico_obs) == "AGE G"] <- "AGE_G"
names(brazil_obs)[names(brazil_obs) == "AGE G"] <- "AGE_G"

mexico_obs <- mexico_obs[mexico_obs$AGE_G == "Total", ]

extract_national_shares <- function(db, state_name, year_value) {
  sub <- db[db$STATE == state_name & db$YEAR == year_value, ]

  c_val <- sum(sub$TOT[sub$RELIGION == "CATOLIC"], na.rm = TRUE)
  p_val <- sum(sub$TOT[sub$RELIGION == "PROTESTANT"], na.rm = TRUE)
  u_val <- sum(sub$TOT[sub$RELIGION == "NONE"], na.rm = TRUE)
  tot <- c_val + p_val + u_val

  c(C = c_val / tot, P = p_val / tot, U = u_val / tot)
}

apply_closed_form <- function(c0, p0, u0, rho_c, rho_p, horizon_years) {
  c0 <- as.numeric(c0)
  p0 <- as.numeric(p0)
  u0 <- as.numeric(u0)
  horizon_decades <- horizon_years / 10

  a_c <- rho_c * horizon_decades
  a_p <- rho_p * horizon_decades
  a_u <- 0

  shift <- max(a_c, a_p, a_u)
  e_c <- exp(a_c - shift)
  e_p <- exp(a_p - shift)
  e_u <- exp(a_u - shift)

  denom <- c0 * e_c + p0 * e_p + u0 * e_u

  c(
    C = as.numeric(c0 * e_c / denom),
    P = as.numeric(p0 * e_p / denom),
    U = as.numeric(u0 * e_u / denom)
  )
}

mexico_anchor_year <- 2020
brazil_anchor_year <- 2022

mexico_anchor <- extract_national_shares(mexico_obs, "Estados Unidos Mexicanos", mexico_anchor_year)
brazil_anchor <- extract_national_shares(brazil_obs, "Brasil", brazil_anchor_year)

base_rows <- db[db[["INTERVAL"]] == "OBS", ]
base_rows <- base_rows[order(base_rows[["ISO3"]]), ]

forecast_rows <- db[db[["INTERVAL"]] != "OBS", ]
forecast_rows <- forecast_rows[order(forecast_rows[["year"]], forecast_rows[["INTERVAL"]], forecast_rows[["ISO3"]]), ]

country_iso3 <- unique(base_rows[["ISO3"]])
country_iso3 <- country_iso3[country_iso3 != "LatAm"]

base_year <- unique(base_rows[["year"]])
if (length(base_year) != 1) {
  stop("The input file must contain a single observed baseline year.")
}
base_year <- base_year[1]

### Fill LB and UB rows for each country.
### Mexico and Brazil use their own national census anchors.
### Other countries keep the 2014 observed baseline with the Mexico/Brazil bounds.
for (iso in country_iso3) {
  base_row <- base_rows[base_rows[["ISO3"]] == iso, ]

  if (nrow(base_row) != 1) {
    stop(paste("Expected one observed baseline row for", iso))
  }

  c0 <- as.numeric(base_row[["C"]])
  p0 <- as.numeric(base_row[["P"]])
  u0 <- as.numeric(base_row[["U"]])

  country_rows <- which(db[["ISO3"]] == iso & db[["INTERVAL"]] %in% c("LB", "UB"))

  for (idx in country_rows) {
    interval_label <- db[idx, "INTERVAL"]

    if (iso == "MEX") {
      shares_now <- apply_closed_form(
        c0 = mexico_anchor["C"],
        p0 = mexico_anchor["P"],
        u0 = mexico_anchor["U"],
        rho_c = rho_mexico_c,
        rho_p = rho_mexico_p,
        horizon_years = as.numeric(db[idx, "year"]) - mexico_anchor_year
      )
    } else if (iso == "BRA") {
      shares_now <- apply_closed_form(
        c0 = brazil_anchor["C"],
        p0 = brazil_anchor["P"],
        u0 = brazil_anchor["U"],
        rho_c = rho_brazil_c,
        rho_p = rho_brazil_p,
        horizon_years = as.numeric(db[idx, "year"]) - brazil_anchor_year
      )
    } else if (interval_label == "LB") {
      shares_now <- apply_closed_form(
        c0 = c0,
        p0 = p0,
        u0 = u0,
        rho_c = rho_mexico_c,
        rho_p = rho_mexico_p,
        horizon_years = as.numeric(db[idx, "year"]) - base_year
      )
    } else {
      shares_now <- apply_closed_form(
        c0 = c0,
        p0 = p0,
        u0 = u0,
        rho_c = rho_brazil_c,
        rho_p = rho_brazil_p,
        horizon_years = as.numeric(db[idx, "year"]) - base_year
      )
    }

    pop_now <- as.numeric(db[idx, "POPULATION"])

    db[idx, "C"] <- shares_now["C"]
    db[idx, "P"] <- shares_now["P"]
    db[idx, "U"] <- shares_now["U"]
    db[idx, "c"] <- shares_now["C"] * pop_now
    db[idx, "p"] <- shares_now["P"] * pop_now
    db[idx, "u"] <- shares_now["U"] * pop_now
  }
}

### Replace the 2014 observed rows for Mexico and Brazil with the model backcast
### from the latest exact census anchor, while keeping the existing file structure.
idx_mex_obs <- which(db[["ISO3"]] == "MEX" & db[["INTERVAL"]] == "OBS")
idx_bra_obs <- which(db[["ISO3"]] == "BRA" & db[["INTERVAL"]] == "OBS")

if (length(idx_mex_obs) != 1 || length(idx_bra_obs) != 1) {
  stop("Expected one observed 2014 row for Mexico and one for Brazil.")
}

mex_obs_shares <- apply_closed_form(
  c0 = mexico_anchor["C"],
  p0 = mexico_anchor["P"],
  u0 = mexico_anchor["U"],
  rho_c = rho_mexico_c,
  rho_p = rho_mexico_p,
  horizon_years = as.numeric(db[idx_mex_obs, "year"]) - mexico_anchor_year
)
bra_obs_shares <- apply_closed_form(
  c0 = brazil_anchor["C"],
  p0 = brazil_anchor["P"],
  u0 = brazil_anchor["U"],
  rho_c = rho_brazil_c,
  rho_p = rho_brazil_p,
  horizon_years = as.numeric(db[idx_bra_obs, "year"]) - brazil_anchor_year
)

db[idx_mex_obs, "C"] <- mex_obs_shares["C"]
db[idx_mex_obs, "P"] <- mex_obs_shares["P"]
db[idx_mex_obs, "U"] <- mex_obs_shares["U"]
db[idx_mex_obs, "c"] <- mex_obs_shares["C"] * as.numeric(db[idx_mex_obs, "POPULATION"])
db[idx_mex_obs, "p"] <- mex_obs_shares["P"] * as.numeric(db[idx_mex_obs, "POPULATION"])
db[idx_mex_obs, "u"] <- mex_obs_shares["U"] * as.numeric(db[idx_mex_obs, "POPULATION"])

db[idx_bra_obs, "C"] <- bra_obs_shares["C"]
db[idx_bra_obs, "P"] <- bra_obs_shares["P"]
db[idx_bra_obs, "U"] <- bra_obs_shares["U"]
db[idx_bra_obs, "c"] <- bra_obs_shares["C"] * as.numeric(db[idx_bra_obs, "POPULATION"])
db[idx_bra_obs, "p"] <- bra_obs_shares["P"] * as.numeric(db[idx_bra_obs, "POPULATION"])
db[idx_bra_obs, "u"] <- bra_obs_shares["U"] * as.numeric(db[idx_bra_obs, "POPULATION"])

### Fill the country-level mean as the simple average between LB and UB.
for (iso in country_iso3) {
  for (yr in sort(unique(db[["year"]][db[["INTERVAL"]] != "OBS"]))) {
    idx_lb <- which(db[["ISO3"]] == iso & db[["year"]] == yr & db[["INTERVAL"]] == "LB")
    idx_ub <- which(db[["ISO3"]] == iso & db[["year"]] == yr & db[["INTERVAL"]] == "UB")
    idx_mean <- which(db[["ISO3"]] == iso & db[["year"]] == yr & db[["INTERVAL"]] == "mean")

    if (length(idx_lb) != 1 || length(idx_ub) != 1 || length(idx_mean) != 1) {
      stop(paste("Expected one LB, UB, and mean row for", iso, "in", yr))
    }

    db[idx_mean, "C"] <- (as.numeric(db[idx_lb, "C"]) + as.numeric(db[idx_ub, "C"])) / 2
    db[idx_mean, "P"] <- (as.numeric(db[idx_lb, "P"]) + as.numeric(db[idx_ub, "P"])) / 2
    db[idx_mean, "U"] <- (as.numeric(db[idx_lb, "U"]) + as.numeric(db[idx_ub, "U"])) / 2
    db[idx_mean, "c"] <- (as.numeric(db[idx_lb, "c"]) + as.numeric(db[idx_ub, "c"])) / 2
    db[idx_mean, "p"] <- (as.numeric(db[idx_lb, "p"]) + as.numeric(db[idx_ub, "p"])) / 2
    db[idx_mean, "u"] <- (as.numeric(db[idx_lb, "u"]) + as.numeric(db[idx_ub, "u"])) / 2
  }
}

### Build the LatAm aggregate from the country forecasts, using population weights.
idx_latam_obs <- which(db[["ISO3"]] == "LatAm" & db[["INTERVAL"]] == "OBS")
idx_country_obs <- which(db[["ISO3"]] %in% country_iso3 & db[["INTERVAL"]] == "OBS")

total_pop_obs <- sum(as.numeric(db[idx_country_obs, "POPULATION"]))
total_c_obs <- sum(as.numeric(db[idx_country_obs, "c"]))
total_p_obs <- sum(as.numeric(db[idx_country_obs, "p"]))
total_u_obs <- sum(as.numeric(db[idx_country_obs, "u"]))

db[idx_latam_obs, "POPULATION"] <- total_pop_obs
db[idx_latam_obs, "c"] <- total_c_obs
db[idx_latam_obs, "p"] <- total_p_obs
db[idx_latam_obs, "u"] <- total_u_obs
db[idx_latam_obs, "C"] <- total_c_obs / total_pop_obs
db[idx_latam_obs, "P"] <- total_p_obs / total_pop_obs
db[idx_latam_obs, "U"] <- total_u_obs / total_pop_obs

forecast_years <- sort(unique(db[["year"]][db[["INTERVAL"]] != "OBS"]))

for (yr in forecast_years) {
  for (interval_label in c("LB", "UB")) {
    idx_latam <- which(db[["ISO3"]] == "LatAm" & db[["year"]] == yr & db[["INTERVAL"]] == interval_label)

    if (length(idx_latam) != 1) {
      stop(paste("Expected one LatAm", interval_label, "row in", yr))
    }

    idx_country <- which(db[["ISO3"]] %in% country_iso3 & db[["year"]] == yr & db[["INTERVAL"]] == interval_label)

    total_pop <- sum(as.numeric(db[idx_country, "POPULATION"]))
    total_c <- sum(as.numeric(db[idx_country, "c"]))
    total_p <- sum(as.numeric(db[idx_country, "p"]))
    total_u <- sum(as.numeric(db[idx_country, "u"]))

    db[idx_latam, "POPULATION"] <- total_pop
    db[idx_latam, "c"] <- total_c
    db[idx_latam, "p"] <- total_p
    db[idx_latam, "u"] <- total_u
    db[idx_latam, "C"] <- total_c / total_pop
    db[idx_latam, "P"] <- total_p / total_pop
    db[idx_latam, "U"] <- total_u / total_pop
  }

  idx_latam_lb <- which(db[["ISO3"]] == "LatAm" & db[["year"]] == yr & db[["INTERVAL"]] == "LB")
  idx_latam_ub <- which(db[["ISO3"]] == "LatAm" & db[["year"]] == yr & db[["INTERVAL"]] == "UB")
  idx_latam_mean <- which(db[["ISO3"]] == "LatAm" & db[["year"]] == yr & db[["INTERVAL"]] == "mean")

  db[idx_latam_mean, "POPULATION"] <- (as.numeric(db[idx_latam_lb, "POPULATION"]) + as.numeric(db[idx_latam_ub, "POPULATION"])) / 2
  db[idx_latam_mean, "c"] <- (as.numeric(db[idx_latam_lb, "c"]) + as.numeric(db[idx_latam_ub, "c"])) / 2
  db[idx_latam_mean, "p"] <- (as.numeric(db[idx_latam_lb, "p"]) + as.numeric(db[idx_latam_ub, "p"])) / 2
  db[idx_latam_mean, "u"] <- (as.numeric(db[idx_latam_lb, "u"]) + as.numeric(db[idx_latam_ub, "u"])) / 2
  db[idx_latam_mean, "C"] <- (as.numeric(db[idx_latam_lb, "C"]) + as.numeric(db[idx_latam_ub, "C"])) / 2
  db[idx_latam_mean, "P"] <- (as.numeric(db[idx_latam_lb, "P"]) + as.numeric(db[idx_latam_ub, "P"])) / 2
  db[idx_latam_mean, "U"] <- (as.numeric(db[idx_latam_lb, "U"]) + as.numeric(db[idx_latam_ub, "U"])) / 2
}

write.csv(db, output_file, row.names = FALSE, fileEncoding = "Windows-1252")
write.csv(db, mirror_output_file, row.names = FALSE, fileEncoding = "Windows-1252")

message("Created: ", output_file)
message("Mirror copy: ", mirror_output_file)
