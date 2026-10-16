#!/usr/bin/env Rscript
# Run in the current container before the first rebuild. The next container
# mounts this host-backed state directory directly; later sessions persist there.
# SQLite's backup API copies live databases consistently, including WAL contents.

suppressPackageStartupMessages({
  library(DBI)
  library(RSQLite)
})

source_dir <- "/root/.codex"
destination_dir <- "/work/.codex_state"
skip_directories <- c(".tmp", ".tmp-private", "tmp", "ipc", "thread-writer-locks", "shell_snapshots", "cache")
dir.create(destination_dir, recursive = TRUE, mode = "0700", showWarnings = FALSE)
stopifnot(normalizePath(source_dir) != normalizePath(destination_dir))
entries <- list.files(source_dir, all.files = TRUE, no.. = TRUE, full.names = TRUE)
for (path in entries) {
  name <- basename(path)
  if (name %in% skip_directories || grepl("[.]sqlite(-wal|-shm)$|[.]lock$", name)) next
  target <- file.path(destination_dir, name)
  if (grepl("[.]sqlite$", name)) {
    source <- dbConnect(SQLite(), path, flags = SQLITE_RO)
    destination <- dbConnect(SQLite(), target)
    tryCatch(sqliteCopyDatabase(source, destination), finally = {
      dbDisconnect(source)
      dbDisconnect(destination)
    })
  } else {
    if (dir.exists(path)) {
      dir.create(target, recursive = TRUE, showWarnings = FALSE)
      files <- list.files(path, recursive = TRUE, all.files = TRUE, no.. = TRUE, full.names = TRUE)
      files <- files[!file.info(files)$isdir]
      for (file in files) {
        relative <- substring(file, nchar(path) + 2L)
        output <- file.path(target, relative)
        dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
        if (!file.copy(file, output, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)) {
          stop("Cannot preserve Codex state file: ", file)
        }
      }
    } else if (!file.copy(path, target, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)) {
      stop("Cannot preserve Codex state file: ", path)
    }
  }
}
message("Codex history and configuration copied to the persistent host directory.")
