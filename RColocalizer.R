#Automated colocalization analysis with RBioformats and BioConductor packages
#PL 2016


#Define function that takes a directory as single argument, analysis will be performed on files inside said directory
RColocalizer <- function(dir) {
  #Set working directory to dir supplied, use directory name as filename for PDF later
  folder_name <- basename(dir)
  setwd(dir)
  
  #Dependencies
  require("rJava")
  require("devtools")
  require("EBImage")
  require("RBioFormats")
  
  #Ask for user input to select channels, assume the same values are held per data set
  chan1 <- as.numeric(readline(prompt = "What is the first channel index for analysis? "))
  fluor1 <- readline(prompt = "What fluorophore was used for the first channel? ")
  chan2 <- as.numeric(readline(prompt = "What is the second channel index for analysis? "))
  fluor2 <- readline(prompt = "What fluorophore was used for the second channel? ")
  
  #Read images
  input_files <- dir(pattern = "cel..oif$")
  images_list <- lapply(input_files, read.image)
  names(images_list) <- paste(input_files, seq_along(images_list))
 
  #Turn images into matrices, this will become a for loop for automated processing of all files
  #Channel 1 = GFP, Channel 2 = mCh, this could become a user input thing once I figure out how that works
  for (i in seq_along(images_list)) {
    
    #Start pdf, plots will go in here
    pdf(paste(folder_name, "_cel", i, ".pdf", sep = ""), paper = "a4")

    chan1_df <- data.frame(getFrame(images_list[[i]], chan1))
    chan2_df <- data.frame(getFrame(images_list[[i]], chan2))
  
    #Change pixel intensity scale from 0 to 1 to 0 to 4096
    if (max(chan1_df, na.rm = TRUE) <= 1) {
      chan1_df_16bit <- chan1_df * 4096
      chan2_df_16bit <- chan2_df * 4096
    }
    
    #Plot histograms to show distributions of intensities
    #Data frames cannot be plotted as histograms, hence as.matrix
    hist(as.matrix(chan1_df_16bit), main = "Pre-normalization", xlab = paste(fluor1, "intensity"))
    hist(as.matrix(chan2_df_16bit), main = "Pre-normalization", xlab = paste(fluor2, "intensity"))
    
    #Scatter plot of un-normalized data
    #Scatter plot takes numeric vectors as input, this line also converts the data frames
    smoothScatter(as.vector(as.matrix(chan1_df_16bit)), as.vector(as.matrix(chan2_df_16bit)), main = "Pre-normalization", xlab = fluor1, ylab = fluor2)
    
    #Convert data frames to vectors for further processing
    chan1_vect <- as.vector(as.matrix(chan1_df_16bit))
    chan2_vect <- as.vector(as.matrix(chan2_df_16bit))
    
    #Exclude values <200 AU
    chan1_vect[chan1_vect < 200] <- NA
    chan2_vect[chan2_vect < 200] <- NA
    
    #Plot new histograms
    hist(chan1_vect, main = "Filter < 200 AU", xlab = paste(fluor1, "intensity"))
    hist(chan2_vect, main = "Filter < 200 AU", xlab = paste(fluor2, "intensity"))
    
    #Normalize channels based on mean fluorescent intensity
    chan1_relative <- chan1_vect / mean(chan1_vect, na.rm = TRUE)
    chan2_relative <- chan2_vect / mean(chan2_vect, na.rm = TRUE)
    
    #Plot histograms, scatter plot, regression line (-1 to force line through origin)
    hist(chan1_relative, main = "Normalized", xlab = paste(fluor1, "intensity"))
    hist(chan2_relative, main = "Normalized", xlab = paste(fluor2, "intensity"))
    smoothScatter(chan1_relative, chan2_relative, main = "Pre-Y axis transformation", xlab = fluor1, ylab = fluor2)
    abline(lm(chan2_relative ~ chan1_relative - 1), col = "red")
    
    #Calculate how Y-axis (mch_trans) needs to be transformed to reach X=Y
    regression_norm <- lm(chan2_relative ~ chan1_relative - 1)
    chan2_trans <- chan2_relative * (1 / regression_norm$coefficients[1])
    
    #Place both vectors in data frame, plot new scatter plot
    scatter_df <- data.frame(chan1_relative, chan2_trans)
    smoothScatter(scatter_df, xlim = c(0, max(scatter_df$chan1_relative, na.rm = TRUE)), ylim = c(0, max(scatter_df$chan2_trans, na.rm = TRUE)), main = "Final", xlab = fluor1, ylab = fluor2)
    abline(lm(chan2_trans ~ chan1_relative - 1, data = scatter_df), col = "red")
    
    #Calculate fineal linear regression
    final_model <- lm(chan2_trans ~ chan1_relative - 1, data = scatter_df)
    dev.off()
    #Create .txt file to output regression results
    sink(file = paste(folder_name, "_cel", i, ".txt", sep = ""), type = "output")
    print(summary(final_model))
    sink()
    #Print to console to keep track of progress
    print(i)
  }
}