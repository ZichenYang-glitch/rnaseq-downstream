#!/usr/bin/env -S Rscript --vanilla

# Bootstrap the locked R library from the reviewed lock artifacts.
#
# This script deliberately has no "install latest" fallback.  The bootstrap
# and primary packages come only from r-sources.lock; every other package is
# restored at the version recorded in renv.lock.

EXPECTED_R_VERSION <- "4.6.1"
EXPECTED_SHA256SUM_VERSION <- "9.5"
EXPECTED_BIOCONDUCTOR_VERSION <- "3.23"
EXPECTED_CONDA_PACKAGES <- c(zlib = "1.3.2", libuv = "1.52.1")
TARGET_LIBRARY_ENV <- "RNASEQ_P0_R_LIBRARY"
TARGET_LIBRARY_MARKER <- ".rnaseq-downstream-p0-library"
BOOTSTRAP_PACKAGES <- c("renv", "BiocManager")
PRIMARY_INSTALL_ORDER <- c(
  "limma", "edgeR", "tximport", "DESeq2", "apeglm", "compcodeR", "airway"
)
EXPECTED_SOURCE_VERSIONS <- c(
  renv = "1.2.4",
  BiocManager = "1.30.27",
  limma = "3.68.0",
  edgeR = "4.10.0",
  tximport = "1.40.0",
  DESeq2 = "1.52.0",
  apeglm = "1.34.0",
  compcodeR = "1.48.0",
  airway = "1.32.0"
)

fail <- function(message) {
  stop(message, call. = FALSE)
}

assert_vanilla <- function() {
  if (!"--vanilla" %in% commandArgs(trailingOnly = FALSE)) {
    fail("Run bootstrap.R with Rscript --vanilla to disable user startup state.")
  }
}

script_path <- function() {
  file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_argument) != 1L) {
    fail("Unable to resolve bootstrap.R from the Rscript command line.")
  }
  normalizePath(sub("^--file=", "", file_argument), mustWork = TRUE)
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

  normalizePath(target, mustWork = FALSE)
}

assert_dedicated_library <- function(target, original_library_paths) {
  protected <- unique(vapply(
    original_library_paths,
    normalizePath,
    FUN.VALUE = character(1),
    mustWork = FALSE
  ))
  if (target %in% protected) {
    fail("Refusing to reconcile an R system or site library.")
  }
  if (identical(target, "/")) {
    fail("Refusing to use the filesystem root as an R library.")
  }
}

prepare_target_library <- function(target) {
  existed <- dir.exists(target)
  if (existed) {
    entries <- list.files(target, all.files = TRUE, no.. = TRUE)
    marker <- file.path(target, TARGET_LIBRARY_MARKER)
    if (length(entries) > 0L && !file.exists(marker)) {
      fail(paste(
        "The target library is nonempty and lacks the P0 ownership marker;",
        "refusing to clean it."
      ))
    }
  } else {
    dir.create(target, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(target) || file.access(target, 2L) != 0L) {
    fail(sprintf("The target library is not writable: %s", target))
  }
  marker <- file.path(target, TARGET_LIBRARY_MARKER)
  if (!file.exists(marker)) {
    writeLines("rnaseq-downstream P0 dedicated R library", marker, useBytes = TRUE)
  }
}

assert_conda_package <- function(prefix, package, version) {
  metadata_directory <- file.path(prefix, "conda-meta")
  if (!dir.exists(metadata_directory)) {
    fail("CONDA_PREFIX does not contain a conda-meta directory.")
  }
  prefix_pattern <- paste0(package, "-", version, "-")
  records <- list.files(metadata_directory, full.names = FALSE)
  matches <- records[
    startsWith(records, prefix_pattern) & endsWith(records, ".json")
  ]
  if (length(matches) != 1L) {
    fail(sprintf("The locked Conda package %s %s is required.", package, version))
  }
}

assert_runtime <- function() {
  actual_r <- paste(R.version$major, R.version$minor, sep = ".")
  if (!identical(actual_r, EXPECTED_R_VERSION)) {
    fail(sprintf("R %s is required; found R %s.", EXPECTED_R_VERSION, actual_r))
  }
  if (!identical(Sys.info()[["sysname"]], "Linux") ||
      !identical(R.version$arch, "x86_64")) {
    fail("The P0 lock is supported only on linux-64 (x86_64 Linux).")
  }

  conda_prefix <- Sys.getenv("CONDA_PREFIX", unset = "")
  if (!nzchar(conda_prefix)) {
    fail("The locked Conda environment must be activated (CONDA_PREFIX is unset).")
  }
  conda_prefix <- normalizePath(conda_prefix, mustWork = TRUE)
  r_home <- normalizePath(R.home(), mustWork = TRUE)
  if (!startsWith(r_home, paste0(conda_prefix, "/"))) {
    fail("R was not launched from the activated locked Conda environment.")
  }
  for (package in names(EXPECTED_CONDA_PACKAGES)) {
    assert_conda_package(conda_prefix, package, EXPECTED_CONDA_PACKAGES[[package]])
  }
  zlib_header <- file.path(conda_prefix, "include", "zlib.h")
  if (!file.exists(zlib_header)) {
    fail("The locked zlib development header is missing from CONDA_PREFIX/include.")
  }
  pkg_config <- Sys.which("pkg-config")
  if (!nzchar(pkg_config) ||
      !startsWith(normalizePath(pkg_config, mustWork = TRUE), paste0(conda_prefix, "/"))) {
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

  sha256sum <- Sys.which("sha256sum")
  if (!nzchar(sha256sum)) {
    fail("sha256sum is required from the locked coreutils package.")
  }
  version_output <- suppressWarnings(system2(
    sha256sum,
    "--version",
    stdout = TRUE,
    stderr = TRUE
  ))
  version_status <- attr(version_output, "status")
  if (is.null(version_status)) {
    version_status <- 0L
  }
  expected_pattern <- sprintf(
    "\\(GNU coreutils\\) %s$",
    gsub("\\.", "\\\\.", EXPECTED_SHA256SUM_VERSION)
  )
  if (version_status != 0L || length(version_output) == 0L ||
      !grepl(expected_pattern, version_output[[1L]])) {
    fail(sprintf(
      "GNU sha256sum %s is required from the locked environment.",
      EXPECTED_SHA256SUM_VERSION
    ))
  }

  sha256sum
}

read_source_manifest <- function(path) {
  manifest <- tryCatch(
    read.delim(
      path,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    error = function(condition) fail(sprintf(
      "Unable to read %s: %s",
      path,
      conditionMessage(condition)
    ))
  )
  expected_columns <- c("package", "version", "repository", "role", "url", "sha256")
  if (!identical(names(manifest), expected_columns)) {
    fail("r-sources.lock has an unexpected schema.")
  }
  expected_names <- c(BOOTSTRAP_PACKAGES, PRIMARY_INSTALL_ORDER)
  if (nrow(manifest) != length(expected_names) ||
      !setequal(manifest$package, expected_names) ||
      anyDuplicated(manifest$package)) {
    fail("r-sources.lock must contain exactly the reviewed source packages.")
  }
  expected_roles <- c(
    renv = "bootstrap",
    BiocManager = "bootstrap",
    limma = "primary",
    edgeR = "primary",
    tximport = "primary",
    DESeq2 = "primary",
    apeglm = "primary",
    compcodeR = "primary",
    airway = "primary"
  )
  actual_roles <- setNames(manifest$role, manifest$package)
  if (!identical(actual_roles[names(expected_roles)], expected_roles)) {
    fail("r-sources.lock contains an unexpected package role.")
  }
  actual_versions <- setNames(manifest$version, manifest$package)
  if (!identical(
    actual_versions[names(EXPECTED_SOURCE_VERSIONS)],
    EXPECTED_SOURCE_VERSIONS
  )) {
    fail("r-sources.lock contains an unexpected package version.")
  }
  expected_repositories <- c(
    renv = "CRAN",
    BiocManager = "CRAN",
    limma = "Bioconductor 3.23",
    edgeR = "Bioconductor 3.23",
    tximport = "Bioconductor 3.23",
    DESeq2 = "Bioconductor 3.23",
    apeglm = "Bioconductor 3.23",
    compcodeR = "Bioconductor 3.23",
    airway = "Bioconductor 3.23"
  )
  actual_repositories <- setNames(manifest$repository, manifest$package)
  if (!identical(
    actual_repositories[names(expected_repositories)],
    expected_repositories
  )) {
    fail("r-sources.lock contains an unexpected package repository.")
  }
  if (any(!grepl("^https://", manifest$url)) ||
      any(!grepl("^[0-9a-f]{64}$", manifest$sha256))) {
    fail("Every source archive must have an HTTPS URL and lowercase SHA-256.")
  }

  manifest
}

verify_archive <- function(path, expected_sha256, sha256sum) {
  output <- suppressWarnings(system2(
    sha256sum,
    c("--", shQuote(path)),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }
  actual_sha256 <- if (length(output) > 0L) {
    strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
  } else {
    ""
  }
  if (status != 0L || !identical(actual_sha256, expected_sha256)) {
    fail(sprintf("SHA-256 verification failed for %s.", basename(path)))
  }
}

file_sha256 <- function(path, sha256sum) {
  output <- suppressWarnings(system2(
    sha256sum,
    c("--", shQuote(path)),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }
  if (status != 0L || length(output) == 0L) {
    fail(sprintf("Unable to hash lock artifact: %s", path))
  }
  digest <- strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
  if (!grepl("^[0-9a-f]{64}$", digest)) {
    fail(sprintf("sha256sum returned an invalid digest for %s.", path))
  }
  digest
}

isolate_renv_state <- function(target, lockfile, manifest_path, sha256sum) {
  fingerprint <- paste0(
    substr(file_sha256(lockfile, sha256sum), 1L, 12L),
    "-",
    substr(file_sha256(manifest_path, sha256sum), 1L, 12L)
  )
  state_root <- paste0(target, "-renv-state-", fingerprint)
  dir.create(state_root, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(state_root) || file.access(state_root, 2L) != 0L) {
    fail(sprintf("The isolated renv state root is not writable: %s", state_root))
  }
  Sys.setenv(RENV_PATHS_ROOT = state_root)
  message(sprintf("Using lock-specific renv state: %s", state_root))
}

download_sources <- function(manifest, destination, sha256sum) {
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  paths <- setNames(character(nrow(manifest)), manifest$package)

  for (index in seq_len(nrow(manifest))) {
    record <- manifest[index, , drop = FALSE]
    archive <- file.path(destination, basename(record$url))
    message(sprintf(
      "Downloading %s %s from the locked source URL.",
      record$package,
      record$version
    ))
    tryCatch(
      download.file(record$url, archive, mode = "wb", quiet = FALSE),
      error = function(condition) fail(sprintf(
        "Download failed for %s %s: %s",
        record$package,
        record$version,
        conditionMessage(condition)
      ))
    )
    verify_archive(archive, record$sha256, sha256sum)
    paths[[record$package]] <- archive
  }

  paths
}

install_verified_source <- function(package, archive, version, target) {
  message(sprintf("Installing verified source archive %s %s.", package, version))
  r_executable <- file.path(R.home("bin"), "R")
  status <- suppressWarnings(system2(
    r_executable,
    c(
      "CMD",
      "INSTALL",
      "--preclean",
      "--no-multiarch",
      paste0("--library=", shQuote(target)),
      shQuote(archive)
    ),
    env = c(
      paste0("R_LIBS=", shQuote(target)),
      paste0("R_LIBS_USER=", shQuote(target)),
      "R_PROFILE_USER=",
      "R_ENVIRON_USER="
    )
  ))
  if (!identical(status, 0L)) {
    fail(sprintf(
      "R CMD INSTALL failed for %s %s with exit status %s.",
      package,
      version,
      status
    ))
  }
  actual <- as.character(packageVersion(package, lib.loc = target))
  if (!identical(actual, version)) {
    fail(sprintf(
      "Installed %s %s, but the manifest requires %s.",
      package,
      actual,
      version
    ))
  }
  loaded <- suppressWarnings(suppressPackageStartupMessages(
    requireNamespace(package, quietly = TRUE, lib.loc = target)
  ))
  if (!isTRUE(loaded)) {
    fail(sprintf("Installed package namespace cannot be loaded: %s.", package))
  }
  namespace_path <- normalizePath(
    getNamespaceInfo(asNamespace(package), "path"),
    mustWork = TRUE
  )
  target_prefix <- paste0(normalizePath(target, mustWork = TRUE), "/")
  if (!startsWith(namespace_path, target_prefix)) {
    fail(sprintf("Installed package was not loaded from the target: %s.", package))
  }
}

run_verifier <- function(root, target) {
  verifier <- file.path(root, "environment", "verify.R")
  rscript <- file.path(R.home("bin"), "Rscript")
  output <- suppressWarnings(system2(
    rscript,
    c("--vanilla", shQuote(verifier), "--library", shQuote(target)),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }
  if (length(output) > 0L) {
    cat(paste(output, collapse = "\n"), "\n", sep = "")
  }
  if (status != 0L) {
    fail("The restored P0 R library failed environment/verify.R.")
  }
}

main <- function() {
  assert_vanilla()
  root <- dirname(dirname(script_path()))
  lockfile <- file.path(root, "renv.lock")
  manifest_path <- file.path(root, "environment", "r-sources.lock")
  if (!file.exists(lockfile) || !file.exists(manifest_path)) {
    fail("renv.lock and environment/r-sources.lock are both required.")
  }

  original_library_paths <- .libPaths()
  target <- parse_target_library(commandArgs(trailingOnly = TRUE))
  assert_dedicated_library(target, original_library_paths)
  sha256sum <- assert_runtime()

  prepare_target_library(target)
  .libPaths(c(target, .Library))
  isolate_renv_state(target, lockfile, manifest_path, sha256sum)

  manifest <- read_source_manifest(manifest_path)
  source_directory <- file.path(tempdir(), "rnaseq-downstream-p0-sources")
  archives <- download_sources(manifest, source_directory, sha256sum)
  versions <- setNames(manifest$version, manifest$package)

  for (package in BOOTSTRAP_PACKAGES) {
    install_verified_source(package, archives[[package]], versions[[package]], target)
  }

  if (!requireNamespace("renv", quietly = TRUE, lib.loc = target)) {
    fail("The verified renv bootstrap package could not be loaded.")
  }
  if (!requireNamespace("BiocManager", quietly = TRUE, lib.loc = target)) {
    fail("The verified BiocManager bootstrap package could not be loaded.")
  }
  options(repos = BiocManager::repositories(version = EXPECTED_BIOCONDUCTOR_VERSION))
  options(renv.config.bioconductor.version = EXPECTED_BIOCONDUCTOR_VERSION)
  message("Refreshing repository availability metadata in the isolated renv state.")
  renv::refresh()

  excluded <- c(BOOTSTRAP_PACKAGES, PRIMARY_INSTALL_ORDER)
  message("Restoring the exact transitive R closure from renv.lock.")
  renv::restore(
    project = root,
    library = target,
    lockfile = lockfile,
    exclude = excluded,
    clean = TRUE,
    strict = TRUE,
    transactional = TRUE,
    retry = FALSE,
    prompt = FALSE
  )

  for (package in PRIMARY_INSTALL_ORDER) {
    install_verified_source(package, archives[[package]], versions[[package]], target)
  }

  run_verifier(root, target)
}

tryCatch(
  main(),
  error = function(condition) {
    message(sprintf("P0 bootstrap failed: %s", conditionMessage(condition)))
    quit(save = "no", status = 1L, runLast = FALSE)
  }
)
