# Function to extract stuffs from options


get.value <- function(NAME, DATA = read.csv2(file.path("./00_Data/00_FileInfos", "Options.csv"))){
  res <- DATA[which(DATA[,1] == NAME), 2]
  res <- as.character(res)  
  return(res)  
  
}
