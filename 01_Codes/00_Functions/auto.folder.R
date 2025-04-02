# Function to create the right list of folder based on Options.csv

# Library

library(tidyverse)

# Function

auto.folder <- function(DATA = read.csv2(file.path("./00_Data/00_FileInfos", "Options.csv"))){
  
  for(x in 1:nrow(DATA)){
    # Vérifier si c'est un path
    if(str_detect(DATA[x,"x"], pattern = "path|log")==TRUE){
      
      PATH <- DATA[x, "value"] %>% as.character()
      
      # Vérifier si le path exist 
      
      if(file.exists(PATH) == F) {
        
        #print(PATH)
        switch(menu(title = paste("Do you want to create this folder?", PATH), graphics = F,
                    choice = c("Yes", "No")
        )+1,
        # Answer 0
        cat("Nothing done ...\n"),
        
        # Answer 1 (yes)
        {
          cat("\nFolder created!\n\n")
          dir.create(file.path(PATH), recursive = T)
          #dir.create(file.path(PATH, "test"), recursive = T) # Don't know why this line exist
        },
        
        # Answer 2 (no)
        cat ("\nFolder not created!\n\n")
        ) # End of the switch menu
        
      }
      
    }
    
  } # End of the loop
  
  cat("\nFolder checking is over!!\n\n")
  
} # End of the function


dir.check <- function(path){
  if(file.exists(path)==F){
    dir.create(path)
    cat(paste("The folder", path, "was created"), sep = "/n")
  } else {
    cat(paste("The folder", path, "is already there, nothing to do"), sep = "/n")
  }
}

