
# Function to Plot the Kegg pathway
plot_pathview <- function(..., save_image = FALSE)
{
  print("eimai sti line5")
  msg <- capture.output(pathview::pathview(...), type = "message")
  print("eimai sti line7")
  msg <- grep("image file", msg, value = T)
  print("eimai sti line8")
  filename <- sapply(strsplit(msg, " "), function(x) x[length(x)])
  print("eimai sti line12")
  img <- png::readPNG(filename)
  print("eimai sti line13")
  grid::grid.raster(img)
  # newFile1=gsub("\\..*", ".png",filename)
  # newFile=gsub("\\..*", ".xml",filename)
  # if(save_image==FALSE) invisible(file.remove(c(filename,newFile1,newFile) ) )
  
}








