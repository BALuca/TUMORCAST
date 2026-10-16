# Install the library used by the analysis image, from its recorded versions.
install.packages("renv", repos = "https://cloud.r-project.org")
renv::restore(
  project = "/opt/tumorcast",
  library = "/opt/tumorcast-library",
  lockfile = "/opt/tumorcast/renv.lock",
  prompt = FALSE
)
