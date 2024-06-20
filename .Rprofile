library(utils)

start_time <- Sys.time()

# Retrieve SLURM job ID and array task ID if they exist
slurm_job_id <- Sys.getenv("SLURM_JOB_ID", unset = "")
slurm_array_task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")

# Define a unique profile file name based on the timestamp and SLURM IDs
profile_file <- paste0("./SlurmLog/renv_startup_profile_",
                       format(start_time, "%Y%m%d_%H%M%S"),
                       "_job", slurm_job_id,
                       "_array", slurm_array_task_id,
                       ".out")

# Start profiling
Rprof(profile_file, line.profiling = TRUE)

# Activate renv
source("renv/activate.R")

# Stop profiling
Rprof(NULL)

end_time <- Sys.time()

# Check if the process took less than 60 seconds
if (as.numeric(difftime(end_time, start_time, units = "secs")) < 60) {
  file.remove(profile_file) # Delete the profile file
}

