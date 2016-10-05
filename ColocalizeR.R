#Automated colocalization analysis with RBioformats and BioConductor packages
#PL 2016


#Define function that takes a directory as single argument, analysis will be performed on files inside said directory
ColocalizeR <- function(dir) {
  #Set working directory to dir supplied, use directory name as filename for PDF later
  folder_name <- basename(dir)
  setwd(dir)
  
  #Dependencies
  require("rJava")
  require("devtools")
  require("EBImage")
  require("RBioFormats")
  require("tools")
  require("data.table")
  require("parallel")
  
  #Ask for user input to select channels, assume the same values are held per data set
  chan1 <- as.numeric(readline(prompt = "What is the first channel index for analysis? "))
  fluor1 <- readline(prompt = "What fluorophore was used for the first channel? ")
  chan2 <- as.numeric(readline(prompt = "What is the second channel index for analysis? "))
  fluor2 <- readline(prompt = "What fluorophore was used for the second channel? ")
  
  #Read image file list
  input_files <- dir(pattern = ".oif$|.lif$|.lsm$|.tif$")
  
  #set up parallel computing for image reading optimization
  total_cores <- detectCores() - 1
  cl <- makeCluster(total_cores)
  
  #set bools to false so they don't break the script
  oif_bool <- FALSE
  lif_bool <- FALSE
  lsm_bool <- FALSE
  tif_bool <- FALSE
  
  if (sum(grepl(".oif$", input_files)) / length(input_files) == 1) {
    images_list <- parLapply(cl, input_files, read.image)
    oif_meta <- list()
    for (i in seq_along(images_list)) {
      oif_meta[[i]] <- globalMetadata(images_list[[i]])
    }
    oif_bool <- TRUE
  } else if (sum(grepl(".lif$", input_files)) / length(input_files) == 1) {
    images_list <- read.image(input_files)
    lif_meta <- list()
    for (i in seq_along(images_list)) {
      lif_meta[[i]] <- seriesMetadata(images_list[[i]])
    }
    lif_bool <- TRUE
  } else if (sum(grepl(".lsm$", input_files)) / length(input_files) == 1) {
    images_list <- parLapply(cl, input_files, read.image)
    lsm_meta <- list()
    for (i in seq_along(images_list)) {
      lsm_meta[[i]] <- seriesMetadata(images_list[[i]])
    }
    lsm_bool <- TRUE
  } else if (sum(grepl(".tif$", input_files)) / length(input_files) == 1) {
    images_list <- parLapply(cl, input_files, read.image)
    tif_meta <- list()
    for (i in seq_along(images_list)) {
      tif_meta[[i]] <- file_path_sans_ext(input_files[[i]])
    }
    tif_bool <- TRUE
  } else {
    stop("Only upload .lif, .oif or .lsm files seperately.")
  }
  

  #Turn images into matrices, below is the function that will be called later with parLapply
  
  colocmainloop <- function(i) {
    require("rJava")
    require("devtools")
    require("EBImage")
    require("RBioFormats")
    require("tools")
    require("data.table")

    if (oif_bool == TRUE) {
      img_name <- file_path_sans_ext(oif_meta[[i]]$`[File Info] DataName`)
    } else if (lif_bool == TRUE) {
      img_name <- lif_meta[[i]]$`Image name`
    } else if (lsm_bool == TRUE) {
      img_name <- lsm_meta[[i]]$`Recording Name #1`
    } else if (tif_bool == TRUE) {
      img_name <- tif_meta[[i]]
    }
    
    #Start pdf, plots will go in here
    pdf(paste(folder_name, "_", i, "_", img_name, ".pdf", sep = ""), paper = "a4")
    
    chan1_dt <- data.table(getFrame(images_list[[i]], chan1))
    chan2_dt <- data.table(getFrame(images_list[[i]], chan2))
    
    #Change pixel intensity scale from 0 to 1 to 0 to 4096
    if (max(chan1_dt, na.rm = TRUE) <= 1) {
      chan1_dt_16bit <- chan1_dt * 4096
      chan2_dt_16bit <- chan2_dt * 4096
    }
    
    #Plot histograms to show distributions of intensities
    #Data frames cannot be plotted as histograms, hence as.matrix
    hist(as.matrix(chan1_dt_16bit), main = "Pre-normalization", xlab = paste(fluor1, "intensity"))
    hist(as.matrix(chan2_dt_16bit), main = "Pre-normalization", xlab = paste(fluor2, "intensity"))
    
    #Scatter plot of un-normalized data
    #Scatter plot takes numeric vectors as input, this line also converts the data frames
    smoothScatter(as.vector(as.matrix(chan1_dt_16bit)), as.vector(as.matrix(chan2_dt_16bit)), main = "Pre-normalization", xlab = fluor1, ylab = fluor2)
    
    #Convert data frames to vectors for further processing
    chan1_vect <- as.vector(as.matrix(chan1_dt_16bit))
    chan2_vect <- as.vector(as.matrix(chan2_dt_16bit))
    
    #Exclude pixels that are < 20% of max in both channels
    chan1_vect[chan1_vect < (0.20 * max(chan1_vect)) & chan2_vect < (0.20 * max(chan2_vect))] <- NA
    
    #Plot new histograms
    hist(chan1_vect, main = "Filter < 20%", xlab = paste(fluor1, "intensity"))
    hist(chan2_vect, main = "Filter < 20%", xlab = paste(fluor2, "intensity"))
    
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
    scatter_dt <- data.table(chan1_relative, chan2_trans)
    smoothScatter(scatter_dt, xlim = c(0, max(scatter_dt$chan1_relative, na.rm = TRUE)), ylim = c(0, max(scatter_dt$chan2_trans, na.rm = TRUE)), main = "Final", xlab = fluor1, ylab = fluor2)
    abline(lm(chan2_trans ~ chan1_relative - 1, data = scatter_dt), col = "red")
    
    #Calculate fineal linear regression
    final_model <- lm(chan2_trans ~ chan1_relative - 1, data = scatter_dt)
    dev.off()
    #Create .txt file to output regression results
    sink(file = paste(folder_name, "_", i, "_", img_name, ".txt", sep = ""), type = "output")
    print(summary(final_model))
    sink()
    #Print to console to keep track of progress
    print(i)
  }
  
  parLapply(cl, seq_along(images_list), colocmainloop)  
  stopCluster(cl)
  
}