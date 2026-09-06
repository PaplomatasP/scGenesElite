# Shared upload parsing, validation, and preview helpers.

SCGENES_MAX_UPLOAD_BYTES <- 100 * 1024^2
SCGENES_MAX_DATA_ROWS <- 50000L
SCGENES_MAX_GENE_COLUMNS <- 50000L
SCGENES_MAX_EXPRESSION_VALUES <- 25000000
SCGENES_MAX_OBJECT_BYTES <- 256 * 1024^2

abort_upload <- function(message) {
  if (requireNamespace("shiny", quietly = TRUE) &&
      !is.null(shiny::getDefaultReactiveDomain())) {
    shiny::validate(shiny::need(FALSE, message))
  }
  stop(message, call. = FALSE)
}

validate_upload_metadata <- function(upload, expected_extension) {
  if (is.null(upload) || is.null(upload$datapath)) {
    abort_upload("Please upload a file before running the analysis.")
  }

  upload_size <- suppressWarnings(as.numeric(upload$size[[1]]))
  if (is.finite(upload_size) && upload_size > SCGENES_MAX_UPLOAD_BYTES) {
    abort_upload(
      sprintf(
        "The uploaded file is %.1f MB. The maximum allowed size is 100 MB.",
        upload_size / 1024^2
      )
    )
  }

  uploaded_extension <- tolower(tools::file_ext(upload$name[[1]]))
  if (!identical(uploaded_extension, tolower(expected_extension))) {
    abort_upload(
      sprintf("Expected a .%s file, but received '%s'.", expected_extension, upload$name[[1]])
    )
  }

  invisible(TRUE)
}

normalize_identifier_column <- function(data) {
  data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(data) < 2L) {
    return(data)
  }

  first_name <- tolower(trimws(colnames(data)[1]))
  identifier_names <- c("", "x", "...1", "row.names", "rownames",
                        "cell", "cell_id", "sample", "sample_id")
  first_values <- data[[1]]
  first_text <- trimws(as.character(first_values))
  looks_like_identifier <-
    first_name %in% identifier_names &&
    !is.numeric(first_values) &&
    !anyNA(first_values) &&
    all(nzchar(first_text)) &&
    !anyDuplicated(first_text)

  if (looks_like_identifier) {
    rownames(data) <- first_text
    data <- data[-1]
  }

  data
}

inspect_expression_dataset <- function(
    data,
    max_rows = SCGENES_MAX_DATA_ROWS,
    max_genes = SCGENES_MAX_GENE_COLUMNS,
    max_expression_values = SCGENES_MAX_EXPRESSION_VALUES,
    max_object_bytes = SCGENES_MAX_OBJECT_BYTES) {
  errors <- character()

  if (!(is.data.frame(data) || is.matrix(data))) {
    return(list(
      valid = FALSE,
      message = "The uploaded object must be a data frame or matrix.",
      data = data
    ))
  }

  data <- normalize_identifier_column(data)
  row_count <- nrow(data)
  column_count <- ncol(data)

  if (row_count < 2L) {
    errors <- c(errors, "At least two cell/sample rows are required.")
  }
  if (column_count < 3L) {
    errors <- c(errors, "At least two gene columns plus one label column are required.")
  }
  if (row_count > max_rows) {
    errors <- c(errors, sprintf("The dataset has %s rows; the limit is %s.", row_count, max_rows))
  }

  if (column_count >= 2L) {
    gene_count <- column_count - 1L
    if (gene_count > max_genes) {
      errors <- c(errors, sprintf("The dataset has %s gene columns; the limit is %s.", gene_count, max_genes))
    }
    if ((as.double(row_count) * as.double(gene_count)) > max_expression_values) {
      errors <- c(
        errors,
        sprintf(
          "The expression matrix has more than %s values and is too large for the public server.",
          format(max_expression_values, big.mark = ",", scientific = FALSE)
        )
      )
    }

    gene_names <- colnames(data)[seq_len(gene_count)]
    trimmed_gene_names <- trimws(gene_names)
    if (anyNA(gene_names) || any(!nzchar(trimmed_gene_names))) {
      errors <- c(errors, "Gene IDs must not be empty or missing.")
    }
    duplicate_ids <- unique(trimmed_gene_names[
      duplicated(tolower(trimmed_gene_names))
    ])
    if (length(duplicate_ids) > 0L) {
      errors <- c(
        errors,
        paste0(
          "Duplicate gene IDs are not allowed: ",
          paste(utils::head(duplicate_ids, 5L), collapse = ", "),
          if (length(duplicate_ids) > 5L) ", ..." else ""
        )
      )
    }
    if (!identical(gene_names, trimmed_gene_names)) {
      errors <- c(errors, "Gene IDs must not contain leading or trailing whitespace.")
    }

    gene_data <- data[seq_len(gene_count)]
    numeric_columns <- vapply(gene_data, is.numeric, logical(1))
    if (!all(numeric_columns)) {
      invalid_names <- names(numeric_columns)[!numeric_columns]
      errors <- c(
        errors,
        paste0(
          "All gene-expression columns must be numeric. Invalid columns: ",
          paste(utils::head(invalid_names, 5L), collapse = ", "),
          if (length(invalid_names) > 5L) ", ..." else ""
        )
      )
    } else {
      non_finite_columns <- vapply(
        gene_data,
        function(column) anyNA(column) || any(!is.finite(column)),
        logical(1)
      )
      if (any(non_finite_columns)) {
        errors <- c(
          errors,
          paste0(
            "Gene-expression values must be finite and non-missing. Invalid columns: ",
            paste(utils::head(names(non_finite_columns)[non_finite_columns], 5L), collapse = ", ")
          )
        )
      }
    }

    labels <- data[[column_count]]
    label_text <- trimws(as.character(labels))
    if (anyNA(labels) || any(!nzchar(label_text))) {
      errors <- c(errors, "The final label column must not contain missing or empty labels.")
    } else {
      label_counts <- table(label_text)
      if (length(label_counts) < 2L) {
        errors <- c(errors, "The final label column must contain at least two classes.")
      }
      if (any(label_counts < 2L)) {
        errors <- c(errors, "Every label class must contain at least two samples.")
      }
    }
  }

  if (as.double(object.size(data)) > max_object_bytes) {
    errors <- c(errors, "The expanded dataset is larger than the 256 MB in-memory limit.")
  }

  list(
    valid = length(errors) == 0L,
    message = paste(unique(errors), collapse = " "),
    data = data
  )
}

validated_expression_dataset <- function(data) {
  inspection <- inspect_expression_dataset(data)
  if (!inspection$valid) {
    abort_upload(inspection$message)
  }
  inspection$data
}

read_uploaded_rds <- function(upload) {
  validate_upload_metadata(upload, "rds")
  data <- tryCatch(
    readRDS(upload$datapath),
    error = function(error) abort_upload(
      paste("The RDS file could not be read:", conditionMessage(error))
    )
  )
  validated_expression_dataset(data)
}

read_uploaded_csv <- function(upload, header = TRUE, sep = ",", quote = "\"") {
  validate_upload_metadata(upload, "csv")
  data <- tryCatch(
    utils::read.table(
      upload$datapath,
      header = isTRUE(as.logical(header)),
      sep = sep,
      quote = quote,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      comment.char = "",
      fill = FALSE
    ),
    error = function(error) abort_upload(
      paste("The CSV file could not be read with the selected options:", conditionMessage(error))
    )
  )
  validated_expression_dataset(data)
}

preview_expression_data <- function(data, max_rows = 10L, max_columns = 11L) {
  row_index <- seq_len(min(nrow(data), max_rows))
  column_index <- tail(seq_len(ncol(data)), min(ncol(data), max_columns))
  data[row_index, column_index, drop = FALSE]
}
