
plot_pathview <- function(..., save_image = FALSE) {
  msg <- capture.output(pathview::pathview(...), type = "message")
  msg <- grep("image file", msg, value = T)
  filename <- sapply(strsplit(msg, " "), function(x) x[length(x)])
  
  cat("filename:", filename, "\n")
  
  # Check if the file exists
  if (!file.exists(filename)) {
    stop("File does not exist:", filename)
  }
  
  # Save the image file to a temporary directory on the server
  tmp_dir <- tempdir()
  img_file <- file.path(tmp_dir, filename)
  file.copy(filename, img_file)
  
  cat("img_file:", img_file, "\n")
  
  # Use the full path to the image file
  img <- png::readPNG(img_file)
  grid::grid.raster(img)
  
  # Optionally remove the image file if save_image is set to FALSE
  if(save_image == FALSE) file.remove(img_file)
}




# 
# plot_pathview <- function(..., save_image = FALSE)
# {
#   msg <- capture.output(pathview::pathview(...), type = "message")
#   msg <- grep("image file", msg, value = T)
#   filename <- sapply(strsplit(msg, " "), function(x) x[length(x)])
#   file.copy(filePath, file)
#   img <- png::readPNG(filename)
#   
#   grid::grid.raster(img)
#   newFile1=gsub("\\..*", ".png",filename)
#   newFile=gsub("\\..*", ".xml",filename)
#   if(save_image==FALSE) invisible(file.remove(c(filename,newFile1,newFile) ) )
#   
# }










