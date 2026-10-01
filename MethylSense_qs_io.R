# MethylSense shared methylation-object reader. Version: 6.0.2. Date: 2026-10-01.
# New files are qs2 (.qs2). Files written by the old qs package (.qs) still open
# when the legacy qs package is installed. The two formats are not interchangeable.

read_methylsense_qs <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("read_methylsense_qs: path must be a single non-empty file path")
  }
  if (!file.exists(path)) {
    stop("read_methylsense_qs: file not found: ", path)
  }

  ext <- tolower(tools::file_ext(path))
  if (ext == "qs2") {
    if (!requireNamespace("qs2", quietly = TRUE)) {
      stop(
        "read_methylsense_qs: package 'qs2' is required to read .qs2 files. ",
        "Install it with install.packages('qs2')."
      )
    }
    return(qs2::qs_read(path))
  }

  if (ext == "qs") {
    if (!requireNamespace("qs", quietly = TRUE)) {
      stop(
        "read_methylsense_qs: package 'qs' is required to read legacy .qs files. ",
        "Install it with install.packages('qs'), or convert the file with qs2::qs_save()."
      )
    }
    return(qs::qread(path))
  }

  stop("read_methylsense_qs: expected a .qs2 or .qs file, got: ", path)
}
