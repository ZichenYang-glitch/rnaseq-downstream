#!/usr/bin/env -S Rscript --vanilla

# Verify the primary P0 runtime contract and emit one compact JSON object.

TARGET_LIBRARY_ENV <- "RNASEQ_P0_R_LIBRARY"
EXPECTED_R_VERSION <- "4.6.1"
EXPECTED_BIOCONDUCTOR_VERSION <- "3.23"
EXPECTED_CONDA_PACKAGES <- c(
    coreutils = "9.5",
    zlib = "1.3.2",
    libuv = "1.52.1",
    numpy = "2.4.4"
)
EXPECTED_PACKAGES <- c(
  renv = "1.2.4",
  BiocManager = "1.30.27",
  BiocVersion = "3.23.1",
  limma = "3.68.0",
  edgeR = "4.10.0",
  tximport = "1.40.0",
  compcodeR = "1.48.0",
  airway = "1.32.0"
)

fail <- function(message) {
  stop(message, call. = FALSE)
}

assert_vanilla <- function() {
  if (!"--vanilla" %in% commandArgs(trailingOnly = FALSE)) {
    fail("Run verify.R with Rscript --vanilla to disable user startup state.")
  }
}

json_string <- function(value) {
  if (anyNA(value)) {
    fail("Cannot encode NA as a JSON string.")
  }
  vapply(as.character(value), function(element) {
    codepoints <- utf8ToInt(enc2utf8(element))
    encoded <- vapply(codepoints, function(codepoint) {
      if (codepoint == 0x22L) {
        "\\\""
      } else if (codepoint == 0x5cL) {
        "\\\\"
      } else if (codepoint < 0x20L) {
        sprintf("\\u%04x", codepoint)
      } else {
        intToUtf8(codepoint)
      }
    }, FUN.VALUE = character(1), USE.NAMES = FALSE)
    paste0("\"", paste0(encoded, collapse = ""), "\"")
  }, FUN.VALUE = character(1), USE.NAMES = FALSE)
}

emit_error <- function(message) {
  document <- paste0(
    '{"schema_version":"1.0","environment":"p0","status":"error",',
    '"errors":[{"code":"ENVIRONMENT_VERIFICATION_FAILED","message":',
    json_string(message),
    "}]}"
  )
  cat(document, "\n", sep = "")
}

parse_target_library <- function(arguments) {
  target <- Sys.getenv(TARGET_LIBRARY_ENV, unset = "")
  index <- 1L

  while (index <= length(arguments)) {
    argument <- arguments[[index]]
    if (identical(argument, "--library")) {
      if (index == length(arguments)) {
        fail("--library requires an absolute path.")
      }
      target <- arguments[[index + 1L]]
      index <- index + 2L
    } else if (startsWith(argument, "--library=")) {
      target <- sub("^--library=", "", argument)
      index <- index + 1L
    } else {
      fail(sprintf("Unknown argument: %s", argument))
    }
  }

  if (!nzchar(target)) {
    fail(sprintf(
      "An explicit target library is required via --library or %s.",
      TARGET_LIBRARY_ENV
    ))
  }
  if (!grepl("^/", target)) {
    fail("The target library must be an absolute path.")
  }
  target <- normalizePath(target, mustWork = FALSE)
  if (!dir.exists(target)) {
    fail(sprintf("The target library does not exist: %s", target))
  }

  target
}

installed_version <- function(package, target) {
  location <- find.package(package, lib.loc = target, quiet = TRUE)
  if (length(location) != 1L || !nzchar(location)) {
    fail(sprintf("Required package is missing from the target library: %s", package))
  }
  as.character(packageVersion(package, lib.loc = target))
}

verify_conda_runtime <- function() {
  prefix <- Sys.getenv("CONDA_PREFIX", unset = "")
  if (!nzchar(prefix)) {
    fail("The locked Conda environment must be activated (CONDA_PREFIX is unset).")
  }
  prefix <- normalizePath(prefix, mustWork = TRUE)
  r_home <- normalizePath(R.home(), mustWork = TRUE)
  if (!startsWith(r_home, paste0(prefix, "/"))) {
    fail("R was not launched from the activated locked Conda environment.")
  }
  metadata_directory <- file.path(prefix, "conda-meta")
  if (!dir.exists(metadata_directory)) {
    fail("CONDA_PREFIX does not contain a conda-meta directory.")
  }
  records <- list.files(metadata_directory, full.names = FALSE)
  for (package in names(EXPECTED_CONDA_PACKAGES)) {
    version <- EXPECTED_CONDA_PACKAGES[[package]]
    prefix_pattern <- paste0(package, "-", version, "-")
    matches <- records[
      startsWith(records, prefix_pattern) & endsWith(records, ".json")
    ]
    if (length(matches) != 1L) {
      fail(sprintf("The locked Conda package %s %s is required.", package, version))
    }
  }
  if (!file.exists(file.path(prefix, "include", "zlib.h"))) {
    fail("The locked zlib development header is missing from CONDA_PREFIX/include.")
  }
  pkg_config <- Sys.which("pkg-config")
  if (!nzchar(pkg_config) ||
      !startsWith(normalizePath(pkg_config, mustWork = TRUE), paste0(prefix, "/"))) {
    fail("pkg-config must come from the activated locked Conda environment.")
  }
  libuv_status <- suppressWarnings(system2(
    pkg_config,
    c("--exists", "libuv"),
    stdout = FALSE,
    stderr = FALSE
  ))
  if (!identical(libuv_status, 0L)) {
    fail("pkg-config could not resolve the locked libuv installation.")
  }
  python <- file.path(prefix, "bin", "python")
  if (!file.exists(python)) {
    fail("The locked Python interpreter is missing from CONDA_PREFIX/bin.")
  }
  numpy_output <- suppressWarnings(system2(
    python,
    c(
      "-I",
      "-c",
      shQuote("import numpy; print(numpy.__version__, end='')")
    ),
    stdout = TRUE,
    stderr = TRUE
  ))
  numpy_status <- attr(numpy_output, "status")
  if (is.null(numpy_status)) {
    numpy_status <- 0L
  }
  if (!identical(numpy_status, 0L) ||
      !identical(paste(numpy_output, collapse = "\n"), "2.4.4")) {
    fail("NumPy 2.4.4 cannot be imported by the locked Python interpreter.")
  }
  EXPECTED_CONDA_PACKAGES
}

load_expected_namespaces <- function(target) {
  for (package in names(EXPECTED_PACKAGES)) {
    loaded <- suppressWarnings(suppressPackageStartupMessages(
      requireNamespace(package, quietly = TRUE, lib.loc = target)
    ))
    if (!isTRUE(loaded)) {
      fail(sprintf("Required package namespace cannot be loaded: %s", package))
    }
    namespace_path <- normalizePath(
      getNamespaceInfo(asNamespace(package), "path"),
      mustWork = TRUE
    )
    target_prefix <- paste0(normalizePath(target, mustWork = TRUE), "/")
    if (!startsWith(namespace_path, target_prefix)) {
      fail(sprintf("Package namespace was not loaded from the target library: %s", package))
    }
  }
}

verify_environment <- function(target) {
  actual_r <- paste(R.version$major, R.version$minor, sep = ".")
  if (!identical(actual_r, EXPECTED_R_VERSION)) {
    fail(sprintf("R %s is required; found R %s.", EXPECTED_R_VERSION, actual_r))
  }
  if (!identical(Sys.info()[["sysname"]], "Linux") ||
      !identical(R.version$arch, "x86_64")) {
    fail("The P0 lock is supported only on linux-64 (x86_64 Linux).")
  }
  actual_conda_packages <- verify_conda_runtime()

  .libPaths(c(target, .Library))
  actual_packages <- vapply(
    names(EXPECTED_PACKAGES),
    installed_version,
    FUN.VALUE = character(1),
    target = target
  )
  mismatched <- names(actual_packages)[actual_packages != EXPECTED_PACKAGES]
  if (length(mismatched) > 0L) {
    details <- vapply(mismatched, function(package) {
      sprintf(
        "%s expected %s, found %s",
        package,
        EXPECTED_PACKAGES[[package]],
        actual_packages[[package]]
      )
    }, FUN.VALUE = character(1))
    fail(paste(details, collapse = "; "))
  }

  load_expected_namespaces(target)

  actual_bioconductor <- as.character(BiocManager::version())
  if (!identical(actual_bioconductor, EXPECTED_BIOCONDUCTOR_VERSION)) {
    fail(sprintf(
      "Bioconductor %s is required; found %s.",
      EXPECTED_BIOCONDUCTOR_VERSION,
      actual_bioconductor
    ))
  }

  package_fields <- paste(
    paste0(json_string(names(actual_packages)), ":", json_string(actual_packages)),
    collapse = ","
  )
  conda_fields <- paste(
    paste0(
      json_string(names(actual_conda_packages)),
      ":",
      json_string(actual_conda_packages)
    ),
    collapse = ","
  )
  paste0(
    '{"schema_version":"1.0","environment":"p0","status":"success",',
    '"platform":"linux-64","library":', json_string(target),
    ',"r_version":', json_string(actual_r),
    ',"bioconductor_version":', json_string(actual_bioconductor),
    ',"conda_packages":{', conda_fields, "}",
    ',"packages":{', package_fields, "}}"
  )
}

status <- 0L
document <- tryCatch(
  {
    assert_vanilla()
    target <- parse_target_library(commandArgs(trailingOnly = TRUE))
    verify_environment(target)
  },
  error = function(condition) {
    status <<- 1L
    emit_error(conditionMessage(condition))
    NULL
  }
)
if (!is.null(document)) {
  cat(document, "\n", sep = "")
}
if (status != 0L) {
  quit(save = "no", status = status, runLast = FALSE)
}
