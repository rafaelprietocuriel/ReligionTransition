setwd("C:/Users/rafae/OneDrive/Documentos/Religion")

suppressPackageStartupMessages({
  library(deSolve)
  library(readr)
})

dir.create("output", showWarnings = FALSE)

decode_u_tags_one <- function(x) {
  if (is.na(x)) {
    return(x)
  }

  hits <- gregexpr("<U\\+[0-9A-Fa-f]{4,6}>", x, perl = TRUE)
  if (length(hits[[1]]) == 1 && hits[[1]][1] == -1) {
    return(x)
  }

  vals <- regmatches(x, hits)[[1]]
  repl <- vapply(
    vals,
    function(tag) {
      code <- sub("^<U\\+([0-9A-Fa-f]{4,6})>$", "\\1", tag, perl = TRUE)
      intToUtf8(strtoi(code, base = 16L))
    },
    character(1)
  )

  regmatches(x, hits) <- list(repl)
  Encoding(x) <- "UTF-8"
  x
}

normalize_text_columns <- function(df) {
  char_cols <- vapply(df, is.character, logical(1))

  if (any(char_cols)) {
    df[char_cols] <- lapply(
      df[char_cols],
      function(col) {
        out <- vapply(col, decode_u_tags_one, character(1), USE.NAMES = FALSE)
        Encoding(out) <- "UTF-8"
        out
      }
    )
  }

  df
}

csv_escape_one <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- enc2utf8(x)
  x <- gsub("\"", "\"\"", x, fixed = TRUE)

  needs_quotes <- grepl("[\",\r\n]", x)
  x[needs_quotes] <- paste0("\"", x[needs_quotes], "\"")
  x
}

write_csv_cp1252 <- function(df, path) {
  df <- normalize_text_columns(df)

  lines <- character(nrow(df) + 1L)
  lines[1] <- paste(csv_escape_one(names(df)), collapse = ",")

  if (nrow(df) > 0) {
    for (i in seq_len(nrow(df))) {
      lines[i + 1L] <- paste(
        csv_escape_one(unlist(df[i, ], use.names = FALSE)),
        collapse = ","
      )
    }
  }

  payload <- paste(lines, collapse = "\r\n")
  payload <- decode_u_tags_one(payload)
  bytes <- iconv(enc2utf8(payload), from = "UTF-8", to = "CP1252", toRaw = TRUE)[[1]]
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(bytes, con)
}

#### country-level parameters with Unaffiliated as the reference group
mexico_rho_C <- -0.183255220747426
mexico_rho_P <-  0.067765122086337

brazil_rho_C <- -0.495485629074094
brazil_rho_P <-  0.0589895036188038


#### raw census inputs
mexico_raw <- as.data.frame(
  readr::read_csv(
    "input/ReligionMexico 20251205.csv",
    show_col_types = FALSE,
    locale = readr::locale(encoding = "UTF-8")
  )
)
brazil_raw <- as.data.frame(
  readr::read_csv(
    "input/ReligionBrazil 20260311.csv",
    show_col_types = FALSE,
    locale = readr::locale(encoding = "UTF-8")
  )
)

mexico_raw <- normalize_text_columns(mexico_raw)
brazil_raw <- normalize_text_columns(brazil_raw)

names(mexico_raw)[1] <- "STATE"
names(brazil_raw)[1] <- "STATE"
names(mexico_raw)[names(mexico_raw) == "AGE G"] <- "AGE_G"
names(brazil_raw)[names(brazil_raw) == "AGE G"] <- "AGE_G"

mexico_raw <- mexico_raw[mexico_raw$AGE_G == "Total", ]


#### simplified state labels used in the output files
simplify_mexico_state <- function(x) {
  x[x == "Estados Unidos Mexicanos"] <- "National"
  x[x == "Coahuila de Zaragoza"] <- "Coahuila"
  x[x == "Michoacán de Ocampo"] <- "Michoacán"
  x[x == "Veracruz de Ignacio de la Llave"] <- "Veracruz"
  x
}

simplify_brazil_state <- function(x) {
  x[x == "Brasil"] <- "National"
  x
}


#### collapse to one row per state-year with exact census shares
build_observed_table <- function(db, country_code, state_mapper) {
  obs <- aggregate(TOT ~ STATE + YEAR + RELIGION, data = db, FUN = sum, na.rm = TRUE)
  obs$STATE <- state_mapper(obs$STATE)
  obs <- obs[order(obs$STATE, obs$YEAR, obs$RELIGION), ]

  states <- unique(obs$STATE)
  out <- vector("list", length(states))

  for (i in seq_along(states)) {
    st <- states[i]
    state_obs <- obs[obs$STATE == st, ]
    years <- sort(unique(state_obs$YEAR))
    rows_i <- vector("list", length(years))

    for (j in seq_along(years)) {
      yr <- years[j]
      sub <- state_obs[state_obs$YEAR == yr, ]

      c_val <- sum(sub$TOT[sub$RELIGION == "CATOLIC"], na.rm = TRUE)
      n_val <- sum(sub$TOT[sub$RELIGION == "NONE"], na.rm = TRUE)
      p_val <- sum(sub$TOT[sub$RELIGION == "PROTESTANT"], na.rm = TRUE)
      pop_val <- c_val + n_val + p_val

      rows_i[[j]] <- data.frame(
        STATE = st,
        YEAR = as.integer(yr),
        C = c_val / pop_val,
        N = n_val / pop_val,
        P = p_val / pop_val,
        COUNTRY = country_code,
        POPULATION = as.integer(round(pop_val)),
        SOURCE = "census",
        stringsAsFactors = FALSE
      )
    }

    out[[i]] <- do.call(rbind, rows_i)
  }

  out <- do.call(rbind, out)
  out[order(out$STATE, out$YEAR), ]
}


#### yearly ODE forecast after the latest census
switching_ode <- function(t, state, parms) {
  C <- state["C"]
  N <- state["N"]
  P <- state["P"]

  rho_C <- parms["rho_C"]
  rho_P <- parms["rho_P"]

  dC <- C * (rho_C * N + (rho_C - rho_P) * P)
  dN <- N * (-rho_C * C - rho_P * P)
  dP <- P * ((rho_P - rho_C) * C + rho_P * N)

  list(c(dC, dN, dP))
}

forecast_state_path <- function(df_state, rho_C, rho_P, end_year = 2050) {
  df_state <- df_state[order(df_state$YEAR), ]
  last_row <- df_state[nrow(df_state), ]
  last_year <- as.integer(last_row$YEAR)

  if (last_year >= end_year) {
    return(data.frame())
  }

  years_out <- seq.int(last_year + 1L, end_year, by = 1L)
  ode_times <- (years_out - last_year) / 10

  sim <- ode(
    y = c(C = last_row$C, N = last_row$N, P = last_row$P),
    times = c(0, ode_times),
    func = switching_ode,
    parms = c(rho_C = rho_C, rho_P = rho_P),
    method = "lsoda"
  )

  sim <- as.data.frame(sim)
  sim <- sim[sim$time > 0, ]

  shares_total <- sim$C + sim$N + sim$P
  sim$C <- pmax(sim$C / shares_total, 0)
  sim$N <- pmax(sim$N / shares_total, 0)
  sim$P <- pmax(sim$P / shares_total, 0)

  data.frame(
    STATE = last_row$STATE,
    YEAR = years_out,
    C = sim$C,
    N = sim$N,
    P = sim$P,
    COUNTRY = last_row$COUNTRY,
    POPULATION = as.integer(round(last_row$POPULATION)),
    SOURCE = "model",
    stringsAsFactors = FALSE
  )
}


#### observed census rows plus yearly model forecast after the latest census
build_country_panel <- function(obs_table, rho_C, rho_P, end_year = 2050) {
  states <- unique(obs_table$STATE)
  out <- vector("list", length(states))

  for (i in seq_along(states)) {
    st <- states[i]
    df_state <- obs_table[obs_table$STATE == st, ]
    fc_state <- forecast_state_path(df_state, rho_C = rho_C, rho_P = rho_P, end_year = end_year)

    if (nrow(fc_state) > 0) {
      out[[i]] <- rbind(df_state, fc_state)
    } else {
      out[[i]] <- df_state
    }
  }

  out <- do.call(rbind, out)
  out[order(out$STATE, out$YEAR), ]
}


#### rebuild the two country panels
mexico_obs <- build_observed_table(mexico_raw, "MEXICO", simplify_mexico_state)
brazil_obs <- build_observed_table(brazil_raw, "BRAZIL", simplify_brazil_state)

mexico_panel <- build_country_panel(
  mexico_obs,
  rho_C = mexico_rho_C,
  rho_P = mexico_rho_P,
  end_year = 2050
)
brazil_panel <- build_country_panel(
  brazil_obs,
  rho_C = brazil_rho_C,
  rho_P = brazil_rho_P,
  end_year = 2050
)

forecast_by_state <- rbind(mexico_panel, brazil_panel)
forecast_by_state <- forecast_by_state[
  order(forecast_by_state$COUNTRY, forecast_by_state$STATE, forecast_by_state$YEAR),
]
row.names(forecast_by_state) <- NULL
forecast_by_state <- normalize_text_columns(forecast_by_state)


#### validation: every census row must match the exact raw input totals
validate_against_input <- function(panel, db, country_code, state_mapper) {
  obs_check <- build_observed_table(db, country_code, state_mapper)
  merged <- merge(
    obs_check,
    panel[panel$SOURCE == "census", ],
    by = c("STATE", "YEAR", "COUNTRY", "SOURCE"),
    suffixes = c("_expected", "_output"),
    all = TRUE
  )

  if (nrow(merged) != nrow(obs_check)) {
    stop(paste("Census row count mismatch for", country_code))
  }

  max_diff <- max(
    abs(merged$C_expected - merged$C_output),
    abs(merged$N_expected - merged$N_output),
    abs(merged$P_expected - merged$P_output),
    na.rm = TRUE
  )

  if (!is.finite(max_diff) || max_diff > 1e-12) {
    stop(paste("Exact census shares do not match for", country_code))
  }
}

validate_against_input(
  forecast_by_state[forecast_by_state$COUNTRY == "MEXICO", ],
  mexico_raw,
  "MEXICO",
  simplify_mexico_state
)
validate_against_input(
  forecast_by_state[forecast_by_state$COUNTRY == "BRAZIL", ],
  brazil_raw,
  "BRAZIL",
  simplify_brazil_state
)

write_csv_cp1252(forecast_by_state, "output/forecast_by_state.csv")

message("Created output/forecast_by_state.csv")
