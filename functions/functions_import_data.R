#Data Import

date_type_column <- function(raw_data){
  
  raw_data$DATE <- gsub(":","/",raw_data$DATE)
  raw_data$DATE <- gsub("Q1","3/31",raw_data$DATE)
  raw_data$DATE <- gsub("Q2","6/30",raw_data$DATE)
  raw_data$DATE <- gsub("Q3","9/30",raw_data$DATE)
  raw_data$DATE <- gsub("Q4","12/31",raw_data$DATE)

  raw_data$DATE <- as.character(raw_data$DATE)
  
  raw_data$DATE <- as.Date(raw_data$DATE, format= "%Y/%m/%d")
  
  return(raw_data)
}

data_import <- function(direction){
  
  sheet = excel_sheets(direction)  
  
  data_1 = read_excel(direction, sheet= sheet[1], na="NaN")
  
  data = lapply(setNames(sheet[2:17], sheet[2:17]),
                function(x) read_excel(direction, sheet=x, na = "NaN")[,2:165])
  
  data = bind_cols(data, .name_repair = "unique")
  
  raw_data = bind_cols(data_1, data)
  
  colnames(raw_data) <- gsub("rinvbf", "RINVBF", colnames(raw_data))
  
  return(date_type_column(raw_data))
}
