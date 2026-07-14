# Generic source-contract validation for publication figures.

value_or <- function(x, y) if (is.null(x)) y else x

contract_states <- function() {
  c(
    "present", "derivable", "compute_required", "decision_required",
    "optional", "blocked"
  )
}

contract_requirements <- function() {
  c("main", "supplement", "optional")
}

contract_privacy_levels <- function() {
  c("public", "internal", "restricted")
}

contract_source_types <- function() {
  c("table", "object", "directory", "none")
}

assert_choice <- function(value, choices, field, contract_id) {
  if (length(value) != 1L || is.na(value) || !value %in% choices) {
    stop(
      sprintf(
        "Contract '%s' has invalid %s '%s'; expected one of: %s",
        contract_id, field, paste(value, collapse = ","),
        paste(choices, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(value)
}

validate_contract_definition <- function(contract) {
  id <- value_or(contract$id, "")
  if (!is.character(id) || length(id) != 1L || !nzchar(id)) {
    stop("Every source contract requires a non-empty string id.", call. = FALSE)
  }

  assert_choice(
    value_or(contract$state, NA_character_),
    contract_states(), "state", id
  )
  assert_choice(
    value_or(contract$required_for, NA_character_),
    contract_requirements(), "required_for", id
  )
  assert_choice(
    value_or(contract$privacy, NA_character_),
    contract_privacy_levels(), "privacy", id
  )

  source <- value_or(contract$source, list(type = "none"))
  assert_choice(
    value_or(source$type, NA_character_),
    contract_source_types(), "source.type", id
  )

  if (source$type != "none") {
    path <- value_or(source$path, "")
    if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
      stop(
        sprintf("Contract '%s' requires source.path.", id),
        call. = FALSE
      )
    }
  }
  if (contract$state == "present" && source$type == "none") {
    stop(
      sprintf("Contract '%s' cannot be present with source.type none.", id),
      call. = FALSE
    )
  }

  columns <- value_or(source$required_columns, character())
  if (!is.character(columns)) {
    stop(
      sprintf("Contract '%s' source.required_columns must be a list.", id),
      call. = FALSE
    )
  }

  invisible(contract)
}

load_contract_registry <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read a contract registry.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("Contract registry not found: ", path, call. = FALSE)
  }

  registry <- yaml::read_yaml(path)
  if (!identical(as.integer(value_or(registry$schema_version, 0L)), 1L)) {
    stop("Contract registry schema_version must be 1.", call. = FALSE)
  }
  contracts <- value_or(registry$contracts, list())
  if (!is.list(contracts) || length(contracts) == 0L) {
    stop("Contract registry must contain a non-empty contracts list.", call. = FALSE)
  }

  invisible(lapply(contracts, validate_contract_definition))
  ids <- vapply(contracts, function(x) x$id, character(1))
  if (anyDuplicated(ids)) {
    stop(
      "Contract ids must be unique: ",
      paste(unique(ids[duplicated(ids)]), collapse = ", "),
      call. = FALSE
    )
  }

  registry
}

resolve_contract_path <- function(project_root, relative_path) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])", relative_path)) {
    stop("Contract paths must be relative to project_root: ", relative_path, call. = FALSE)
  }

  root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  target <- normalizePath(
    file.path(root, relative_path),
    winslash = "/",
    mustWork = FALSE
  )
  if (!identical(target, root) && !startsWith(target, paste0(root, "/"))) {
    stop("Contract path escapes project_root: ", relative_path, call. = FALSE)
  }
  target
}

read_table_columns <- function(path, format = NULL) {
  format <- tolower(value_or(format, tools::file_ext(path)))
  reader <- switch(
    format,
    csv = function() utils::read.csv(
      path, nrows = 0L, check.names = FALSE, stringsAsFactors = FALSE
    ),
    tsv = function() utils::read.delim(
      path, nrows = 0L, check.names = FALSE, stringsAsFactors = FALSE
    ),
    txt = function() utils::read.delim(
      path, nrows = 0L, check.names = FALSE, stringsAsFactors = FALSE
    ),
    stop("Unsupported table format: ", format, call. = FALSE)
  )
  names(reader())
}

table_has_data_row <- function(path, format = NULL) {
  format <- tolower(value_or(format, tools::file_ext(path)))
  reader <- switch(
    format,
    csv = function() utils::read.csv(
      path, nrows = 1L, check.names = FALSE, stringsAsFactors = FALSE
    ),
    tsv = function() utils::read.delim(
      path, nrows = 1L, check.names = FALSE, stringsAsFactors = FALSE
    ),
    txt = function() utils::read.delim(
      path, nrows = 1L, check.names = FALSE, stringsAsFactors = FALSE
    ),
    stop("Unsupported table format: ", format, call. = FALSE)
  )
  nrow(reader()) > 0L
}

audit_source_contract <- function(contract, project_root = ".") {
  validate_contract_definition(contract)

  id <- contract$id
  state <- contract$state
  required_for <- contract$required_for
  source <- value_or(contract$source, list(type = "none"))
  source_type <- source$type
  relative_path <- value_or(source$path, "")
  resolved_path <- ""
  exists <- FALSE
  nonempty <- FALSE
  columns_ok <- NA
  missing_columns <- character()
  error_message <- ""

  if (source_type != "none") {
    resolved_path <- tryCatch(
      resolve_contract_path(project_root, relative_path),
      error = function(e) {
        error_message <<- conditionMessage(e)
        ""
      }
    )
  }

  if (nzchar(resolved_path)) {
    exists <- if (source_type == "directory") {
      dir.exists(resolved_path)
    } else {
      file.exists(resolved_path)
    }
    nonempty <- if (!exists) {
      FALSE
    } else if (source_type == "directory") {
      length(list.files(resolved_path, all.files = FALSE)) > 0L
    } else if (source_type == "table") {
      tryCatch(
        table_has_data_row(resolved_path, source$format),
        error = function(e) {
          error_message <<- conditionMessage(e)
          FALSE
        }
      )
    } else {
      isTRUE(file.info(resolved_path)$size > 0)
    }

    required_columns <- value_or(source$required_columns, character())
    if (source_type == "table" && exists && length(required_columns) > 0L) {
      observed_columns <- tryCatch(
        read_table_columns(resolved_path, source$format),
        error = function(e) {
          error_message <<- conditionMessage(e)
          character()
        }
      )
      missing_columns <- setdiff(required_columns, observed_columns)
      columns_ok <- length(missing_columns) == 0L
    } else if (source_type == "table" && exists) {
      columns_ok <- TRUE
    }
  }

  source_valid <- source_type == "none" ||
    (exists && nonempty && !identical(columns_ok, FALSE) && !nzchar(error_message))
  ready <- identical(state, "present") && source_valid
  blocking <- required_for == "main" && !ready && state != "optional"

  reason <- if (ready) {
    "ready"
  } else if (nzchar(error_message)) {
    error_message
  } else if (state == "present" && !exists) {
    "declared present but source is missing"
  } else if (state == "present" && !nonempty) {
    "declared present but source is empty"
  } else if (state == "present" && identical(columns_ok, FALSE)) {
    paste0("missing columns: ", paste(missing_columns, collapse = ", "))
  } else {
    paste0("state requires action: ", state)
  }

  data.frame(
    id = id,
    state = state,
    required_for = required_for,
    privacy = contract$privacy,
    source_type = source_type,
    source_path = relative_path,
    exists = exists,
    nonempty = nonempty,
    columns_ok = if (is.na(columns_ok)) "" else as.character(columns_ok),
    missing_columns = paste(missing_columns, collapse = ";"),
    ready = ready,
    blocking = blocking,
    producer = value_or(contract$producer, ""),
    reason = reason,
    stringsAsFactors = FALSE
  )
}

audit_contract_registry <- function(registry, project_root = ".") {
  rows <- lapply(
    registry$contracts,
    audit_source_contract,
    project_root = project_root
  )
  do.call(rbind, rows)
}

write_contract_audit <- function(audit, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(audit, path, row.names = FALSE, na = "")
  invisible(path)
}
