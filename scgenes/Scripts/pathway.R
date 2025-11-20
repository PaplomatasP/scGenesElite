
plot_pathview <- function(..., save_image = FALSE) {
  
  # GC Petros S 1: Get working directory before calling pathview
  current_wd <- getwd()
  
  # GC Petros S 2: Extract pathway.id and species from arguments to construct filename
  args_list <- list(...)
  pathway_id <- args_list$pathway.id
  species <- args_list$species
  out_suffix <- ifelse(is.null(args_list$out.suffix), "", args_list$out.suffix)
  
  # GC Petros S 3: Construct expected filename pattern
  # pathview typically creates files like: species_pathway_id.png and species_pathway_id.xml
  if (!is.null(pathway_id) && !is.null(species)) {
    expected_base <- paste0(species, pathway_id)
    expected_png <- paste0(expected_base, ifelse(out_suffix == "", "", paste0(".", out_suffix)), ".pathview.png")
    expected_xml <- paste0(expected_base, ifelse(out_suffix == "", "", paste0(".", out_suffix)), ".xml")
  } else {
    expected_png <- NULL
    expected_xml <- NULL
  }
  
  # GC Petros S 4: Run pathview and capture messages
  # Get list of PNG files before running pathview
  png_files_before <- list.files(pattern = "\\.pathview\\.png$", full.names = FALSE)
  
  # Run pathview and capture any messages (for debugging)
  msg <- suppressWarnings(capture.output({
    pathview::pathview(...)
  }, type = "message"))
  
  # GC Petros S 5: Try to extract filename from messages
  msg_filtered <- grep("image file|png file|Pathview|Working in directory", msg, value = TRUE, ignore.case = TRUE)
  
  filename <- NULL
  if (length(msg_filtered) > 0) {
    # Try multiple patterns to extract filename
    for (msg_line in msg_filtered) {
      # Pattern 1: "image file ... path/to/file.png"
      if (grepl("image file|png file", msg_line, ignore.case = TRUE)) {
        # Try to extract file path
        # Pattern: "image file saved in .../filename.png"
        png_match <- regmatches(msg_line, regexpr("[^\\s]+\\.pathview\\.png", msg_line))
        if (length(png_match) > 0) {
          filename <- png_match[1]
          # Clean up filename (remove quotes, commas, etc.)
          filename <- gsub("['\",;]", "", filename)
          filename <- basename(filename)  # Just use filename, not full path
          break
        }
        
        # Alternative pattern: extract from spaces
        parts <- strsplit(msg_line, "\\s+")[[1]]
        png_idx <- grep("\\.pathview\\.png", parts)
        if (length(png_idx) > 0) {
          filename <- basename(parts[png_idx[1]])
          filename <- gsub("['\",;]", "", filename)
          break
        }
      }
    }
  }
  
  # GC Petros S 5.5: Get list of PNG files after running pathview
  png_files_after <- list.files(pattern = "\\.pathview\\.png$", full.names = FALSE)
  
  # Find newly created files
  new_png_files <- setdiff(png_files_after, png_files_before)
  
  # GC Petros S 6: If filename not found in messages, use newly created files
  if ((is.null(filename) || length(filename) == 0 || filename == "")) {
    # First, try newly created files
    if (length(new_png_files) > 0) {
      filename <- new_png_files[1]  # Use first new file
    } else if (!is.null(expected_png) && file.exists(expected_png)) {
      # Try expected pattern
      filename <- expected_png
    } else {
      # List all .pathview.png files in current directory
      all_pathview_png <- list.files(pattern = "\\.pathview\\.png$", full.names = FALSE)
      
      if (length(all_pathview_png) > 0) {
        # Find most recently created .pathview.png file
        png_info <- file.info(all_pathview_png)
        most_recent_idx <- which.max(png_info$mtime)
        if (length(most_recent_idx) > 0) {
          filename <- all_pathview_png[most_recent_idx]
        }
      }
    }
  }
  
  # GC Petros S 7: Check if filename was found
  if (is.null(filename) || length(filename) == 0 || filename == "" || !file.exists(filename)) {
    # Try alternative: look for any .pathview.png files
    pathview_files <- list.files(pattern = "\\.pathview\\.png$", full.names = FALSE)
    if (length(pathview_files) > 0) {
      # Get most recent
      if (length(pathview_files) > 1) {
        file_info <- file.info(pathview_files)
        filename <- rownames(file_info)[which.max(file_info$mtime)]
      } else {
        filename <- pathview_files[1]
      }
    }
  }
  
  # GC Petros S 8: Final check and error handling
  if (is.null(filename) || length(filename) == 0 || filename == "" || !file.exists(filename)) {
    # List all PNG files for debugging
    all_files <- list.files(pattern = "\\.png$", full.names = FALSE)
    
    if (length(all_files) == 0) {
      stop("Pathview output file not found. The pathway visualization may have failed. ",
           "Please check: 1) Pathway ID is correct, 2) Species matches your data, ",
           "3) Gene IDs match KEGG format (e.g., Entrez IDs).")
    }
    
    # GC Petros S 8.5: Smart file selection - prefer files matching expected pattern
    # Filter files that match expected pattern or are .pathview.png files
    if (!is.null(expected_png)) {
      expected_basename <- basename(expected_png)
      matching_files <- all_files[grepl(gsub("\\.", "\\\\.", expected_basename), all_files, ignore.case = TRUE)]
      
      if (length(matching_files) > 0) {
        # Prefer .pathview.png files
        pathview_files <- matching_files[grepl("\\.pathview\\.png$", matching_files, ignore.case = TRUE)]
        if (length(pathview_files) > 0) {
          filename <- pathview_files[1]
        } else {
          # Prefer files without double dots or extra characters
          clean_files <- matching_files[!grepl("\\.\\..*\\.", matching_files)]
          if (length(clean_files) > 0) {
            filename <- clean_files[1]
          } else {
            filename <- matching_files[1]
          }
        }
      } else {
        # No matching pattern, prefer .pathview.png files
        pathview_files <- all_files[grepl("\\.pathview\\.png$", all_files, ignore.case = TRUE)]
        if (length(pathview_files) > 0) {
          # Get most recent .pathview.png file
          if (length(pathview_files) > 1) {
            file_info <- file.info(pathview_files)
            filename <- rownames(file_info)[which.max(file_info$mtime)]
          } else {
            filename <- pathview_files[1]
          }
        } else {
          # Prefer files without double dots
          clean_files <- all_files[!grepl("\\.\\..*\\.", all_files)]
          if (length(clean_files) > 0) {
            filename <- clean_files[1]
          } else {
            filename <- all_files[1]
          }
        }
      }
    } else {
      # No expected pattern, prefer .pathview.png files
      pathview_files <- all_files[grepl("\\.pathview\\.png$", all_files, ignore.case = TRUE)]
      if (length(pathview_files) > 0) {
        filename <- pathview_files[1]
      } else {
        # Prefer files without double dots (clean filenames)
        clean_files <- all_files[!grepl("\\.\\..*\\.", all_files)]
        if (length(clean_files) > 0) {
          filename <- clean_files[1]
        } else {
          filename <- all_files[1]
        }
      }
    }
    
    # Final check - if still no file found, stop with error
    if (is.null(filename) || length(filename) == 0 || filename == "" || !file.exists(filename)) {
      stop("The filename is not found. Pathview may not have created the output file.")
    }
    
    # Only print info message (not warning) if file was found but not auto-detected
    # This is normal - pathview may create files with slightly different names
  }
  
  # GC Petros S 8.7: Ensure we have a clean filename (remove double dots if any)
  if (grepl("\\.\\..*\\.", filename)) {
    # Try to clean up the filename - remove double dots
    clean_filename <- gsub("\\.{2,}", ".", filename)
    if (file.exists(clean_filename)) {
      filename <- clean_filename
    } else {
      # If cleaning didn't work, try to find a better match
      all_png <- list.files(pattern = paste0("^", gsub("\\.png.*", "", filename), ".*\\.png$"), full.names = FALSE)
      clean_png <- all_png[!grepl("\\.\\..*\\.", all_png)]
      if (length(clean_png) > 0) {
        filename <- clean_png[1]
      }
    }
  }
  
  print(paste("Using pathview output file:", filename))
  
  # GC Petros S 9: Copy to temporary directory and display
  tmp_dir <- tempdir()
  img_file <- file.path(tmp_dir, basename(filename))
  
  # Copy file (overwrite if exists)
  if (file.exists(filename)) {
    file.copy(filename, img_file, overwrite = TRUE)
    
    # GC Petros S 10: Read and display image
    if (file.exists(img_file)) {
      img <- png::readPNG(img_file)
      grid::grid.raster(img)
    } else {
      stop("Failed to copy image file to temporary directory.")
    }
  } else {
    stop(paste("Pathview output file does not exist:", filename))
  }
  
  # GC Petros S 11: Clean up files if save_image is FALSE
  if (save_image == FALSE) {
    # Clean up original files
    newFile1 <- gsub("\\..*", ".png", filename)
    newFile <- gsub("\\..*", ".xml", filename)
    
    # Also try .pathview.png pattern
    if (!file.exists(newFile1) && file.exists(filename)) {
      newFile1 <- filename
    }
    
    # Remove files silently (may not exist)
    files_to_remove <- c(filename, newFile1, newFile)
    files_to_remove <- files_to_remove[file.exists(files_to_remove)]
    
    if (length(files_to_remove) > 0) {
      tryCatch({
        invisible(file.remove(files_to_remove))
      }, error = function(e) {
        # Silently ignore cleanup errors
        invisible(NULL)
      })
    }
  }
}


