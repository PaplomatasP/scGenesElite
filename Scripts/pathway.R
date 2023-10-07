
plot_pathview <- function(..., save_image = FALSE) {


  msg <- capture.output(pathview::pathview(...), type = "message")
  msg <- grep("image file", msg, value = T)
  filename <- sapply(strsplit(msg, " "), function(x) x[length(x)])
  print(filename)
  # Check if the file exists
  full_file_path <- file.path(getwd(), filename)
  if (length(filename) == 0) {
    stop("The filename is not found.")
  }

  # Save the image file to a temporary directory on the server
  tmp_dir <- tempdir()
  img_file <- file.path(tmp_dir, filename)
  file.copy(filename, img_file)

  # Use the full path to the image file
  img <- png::readPNG(img_file)

  grid::grid.raster(img)
  newFile1=gsub("\\..*", ".png",filename)
  newFile=gsub("\\..*", ".xml",filename)
  # Optionally remove the image file if save_image is set to FALSE
  if(save_image==FALSE) invisible(file.remove(c(filename,newFile1,newFile) ) )
}


