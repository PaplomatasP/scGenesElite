
plot_pathview <- function(..., save_image = FALSE) {
  msg <- capture.output(pathview::pathview(...), type = "message")
  msg <- grep("image file", msg, value = T)
  filename <- sapply(strsplit(msg, " "), function(x) x[length(x)])
  
  # Check if the file exists
  if (!file.exists(filename)) {
    stop("File does not exist:", filename)
  }
  
  # Save the image file to a temporary directory on the server
  tmp_dir <- tempdir()
  img_file <- file.path(tmp_dir, filename)
  file.copy(filename, img_file)
  
  # Use the full path to the image file
  img <- png::readPNG(img_file)
  grid::grid.raster(img)
  
  # Optionally remove the image file if save_image is set to FALSE
  if(save_image == FALSE) file.remove(img_file)
}



# # 
# plot_pathview <- function(..., save_image = FALSE)
# {
#   print("eimai sti line5")
#   msg <<- capture.output(pathview::pathview(...), type = "message")
#   print("eimai sti line7")
#   msg <- grep("image file", msg, value = T)
#   print("eimai sti line8")
#   filename <- sapply(strsplit(msg, " "), function(x) x[length(x)])
#   print("eimai sti line12")
#   img <- png::readPNG(filename)
#   print("eimai sti line13")
#   grid::grid.raster(img)
#   newFile1=gsub("\\..*", ".png",filename)
#   newFile=gsub("\\..*", ".xml",filename)
#   if(save_image==FALSE) invisible(file.remove(c(filename,newFile1,newFile) ) )
# 
# }









