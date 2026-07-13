#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
strict <- "--strict" %in% args
args <- args[args != "--strict"]

if (length(args) < 1L || length(args) > 3L) {
  stop(
    paste(
      "Usage: Rscript audit_sources.R <registry.yaml>",
      "[project_root] [output.csv] [--strict]"
    ),
    call. = FALSE
  )
}

command <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command, value = TRUE)
script_path <- normalizePath(
  sub("^--file=", "", file_arg[1]),
  winslash = "/",
  mustWork = TRUE
)
module_root <- dirname(script_path)

source(file.path(module_root, "R", "contracts.R"))

registry_path <- args[1]
project_root <- if (length(args) >= 2L) args[2] else getwd()
output_path <- if (length(args) >= 3L) {
  args[3]
} else {
  file.path(project_root, "results", "publish_figures", "contract_audit.csv")
}

registry <- load_contract_registry(registry_path)
audit <- audit_contract_registry(registry, project_root)
write_contract_audit(audit, output_path)

print(
  audit[, c("id", "state", "required_for", "ready", "blocking", "reason")],
  row.names = FALSE
)
cat(sprintf("\nContract audit written to: %s\n", output_path))

if (strict && any(audit$blocking)) {
  quit(status = 1L)
}
