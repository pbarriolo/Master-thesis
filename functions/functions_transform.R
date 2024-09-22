# functions to select and transform macro data

select_data <- function(df, tcode, year_quarter, forecast_horizon) {
  subset_cols <- c("DATE")

  true_values_index <- forecast_horizon + 1 + which(colnames(df) == paste("CPI", year_quarter, sep = ""))

  for (col in colnames(df)) {
    if (grepl(year_quarter, col)) {
      subset_cols <- append(subset_cols, col)
    }
  }

  df_2 <- df[, subset_cols]


  colnames(df_2) <- c("DATE", colnames(tcode))

  df_2["CPI_1"] <- lag(df_2[["CPI"]], 1)

  df_2[paste("CPI_h", as.character(forecast_horizon), sep = "")] <- lead(df_2[["CPI"]], forecast_horizon + 1)

  df_2["true_CPI"] <- df[, true_values_index]


  return(df_2)
}

transform_column <- function(df, variable, tcode) {
  #### This function stationarizes all series, following the transformation code present in the df tcode
  ####

  df_2 <- df

  if (grepl("CPI", variable) == TRUE) {
    df_2[[variable]] <- 400 * (log(df[[variable]]) - lag(log(df[[variable]]), 1)) # transformation en log diff
  } else if (tcode[[variable]] == 4) {
    df_2[[variable]] <- log(df[[variable]]) # transformation en log
  } else if (tcode[[variable]] == 5) {
    df_2[[variable]] <- 400 * (log(df[[variable]]) - lag(log(df[[variable]]), 1)) # transformation en log diff
  } else if (tcode[[variable]] == 2) {
    df_2[[variable]] <- df[[variable]] - lag(df[[variable]], 1) # transformation en log diff
  }

  return(df_2[variable])
}

transform_data <- function(df, tcode) {
  trans_df <- lapply(colnames(df)[2:length(df)], function(x) transform_column(df, x, tcode))

  temp_df <- do.call(bind_cols, trans_df)

  new_df <- bind_cols(df["DATE"], temp_df)

  colnames(new_df) <- colnames(df)

  # new_df["CPI_1"] <- lag(new_df[["CPI"]], 1)


  return(new_df)
}

get_data <- function(vintage_date, horizon) {
  raw_data <- data_import("C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/data/ALL_DATA.xlsx")

  tcode <- data.frame(CPI = 5, IP = 5, HSTARTS = 4, CUM = 5, M1 = 5, RCOND = 5, RCONS = 5, RG = 5, RINVBF = 5, ROUTPUT = 5, RUC = 2, ULC = 5, WSD = 5, DYS = 1, NAPM = 1, NAPMII = 1, NAPMNOI = 1)

  subset_raw_data <- select_data(raw_data, tcode, vintage_date, horizon)

  # transform series
  df <- transform_data(subset_raw_data[2:nrow(subset_raw_data), ], tcode)

  return(df)
}