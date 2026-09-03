"""Strict request-time evidence capture for local frozen gene-set sources.

This module validates only the public request schema and the identity of the
two declared source byte streams.  GMT and annotation contents are deliberately
not interpreted here; semantic parsing belongs to the later materialization
boundary that copies the verified sources into the private backend workspace.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
import re
import stat
from typing import Any, Mapping

from .errors import InputIntegrityError, InputReadError, InvalidRequestError
from .provenance import require_expected_keys, require_nonempty_string


GENE_SET_IDENTIFIER_TYPE = "symbol"
GENE_SET_GENE_ID_COLUMN = "gene_id"
GENE_SET_SYMBOL_COLUMN = "symbol"

_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
_URI_PATTERN = re.compile(r"^[A-Za-z][A-Za-z0-9+.-]*://")
_GENE_SETS_KEYS = {"gmt", "annotation", "minimum_tested_genes"}
_GMT_KEYS = {"path", "sha256", "collection", "version", "identifier_type"}
_ANNOTATION_KEYS = {
    "path",
    "sha256",
    "name",
    "version",
    "gene_id_column",
    "symbol_column",
}


def _canonical_sha256(value: Any, *, field: str) -> str:
    if not isinstance(value, str) or _SHA256_PATTERN.fullmatch(value) is None:
        raise InvalidRequestError(
            f"'{field}' must be a canonical lowercase SHA-256 digest.",
            details={"field": field, "observed": value},
        )
    return value


def _capture_source(
    declared_path: Any,
    *,
    expected_sha256: Any,
    request_path: Path,
    role: str,
) -> dict[str, Any]:
    """Capture one non-symlink regular UTF-8 file from a stable byte stream."""

    declared = require_nonempty_string(declared_path, field=f"gene_sets.{role}.path")
    if _URI_PATTERN.match(declared):
        raise InvalidRequestError(
            "Gene-set sources must be local files; network URLs are not permitted.",
            details={"field": f"gene_sets.{role}.path", "declared_path": declared},
        )
    digest_expected = _canonical_sha256(
        expected_sha256, field=f"gene_sets.{role}.sha256"
    )
    candidate = Path(declared)
    if not candidate.is_absolute():
        candidate = request_path.parent / candidate

    try:
        lexical = candidate.lstat()
        if stat.S_ISLNK(lexical.st_mode):
            raise OSError("symbolic links are not accepted as gene-set sources")
        if not stat.S_ISREG(lexical.st_mode):
            raise OSError("expected a regular file")
        resolved = candidate.resolve(strict=True)
        before = resolved.stat()
        if not stat.S_ISREG(before.st_mode):
            raise OSError("expected a regular file")
        with resolved.open("rb") as handle:
            content = handle.read()
        after = resolved.stat()
        lexical_after = candidate.lstat()
    except OSError as error:
        raise InputReadError(
            "A declared gene-set source could not be captured as a regular file.",
            path=candidate,
            operation="capture_gene_set_source",
            cause=error,
            details={"role": role, "declared_path": declared},
        ) from error
    except RuntimeError as error:
        raise InputReadError(
            "A declared gene-set source path could not be resolved.",
            path=candidate,
            operation="resolve_gene_set_source",
            cause=error,
            details={"role": role, "declared_path": declared},
        ) from error

    identity_before = (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
    )
    identity_after = (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    )
    lexical_identity_before = (
        lexical.st_dev,
        lexical.st_ino,
        lexical.st_mode,
        lexical.st_size,
        lexical.st_mtime_ns,
    )
    lexical_identity_after = (
        lexical_after.st_dev,
        lexical_after.st_ino,
        lexical_after.st_mode,
        lexical_after.st_size,
        lexical_after.st_mtime_ns,
    )
    if (
        identity_before != identity_after
        or lexical_identity_before != lexical_identity_after
        or stat.S_ISLNK(lexical_after.st_mode)
        or len(content) != after.st_size
    ):
        raise InputIntegrityError(
            "A declared gene-set source changed while it was captured.",
            details={
                "role": role,
                "path": str(resolved),
                "source_identity_stable": identity_before == identity_after,
                "declared_path_identity_stable": lexical_identity_before
                == lexical_identity_after,
                "stat_size_bytes": after.st_size,
                "captured_size_bytes": len(content),
            },
        )

    try:
        content.decode("utf-8")
    except UnicodeDecodeError as error:
        raise InputReadError(
            "A declared gene-set source is not valid UTF-8.",
            path=resolved,
            operation="decode_gene_set_source",
            cause=error,
            details={"role": role, "declared_path": declared},
        ) from error

    digest_observed = hashlib.sha256(content).hexdigest()
    if digest_observed != digest_expected:
        raise InputIntegrityError(
            "A gene-set source does not match its declared SHA-256 digest.",
            details={
                "role": role,
                "path": str(resolved),
                "expected_sha256": digest_expected,
                "observed_sha256": digest_observed,
            },
        )
    return {
        "path": str(resolved),
        "declared_path": declared,
        "sha256": digest_observed,
        "size_bytes": len(content),
    }


def parse_gene_sets_request(
    value: Any,
    *,
    request_path: Path,
) -> dict[str, Any]:
    """Validate and normalize the optional public ``gene_sets`` request object."""

    if not isinstance(value, Mapping):
        raise InvalidRequestError("'gene_sets' must be an object.")
    require_expected_keys(
        value,
        allowed=_GENE_SETS_KEYS,
        required=_GENE_SETS_KEYS,
        context="analysis gene_sets request",
    )

    gmt = value["gmt"]
    if not isinstance(gmt, Mapping):
        raise InvalidRequestError("'gene_sets.gmt' must be an object.")
    require_expected_keys(
        gmt,
        allowed=_GMT_KEYS,
        required=_GMT_KEYS,
        context="analysis gene_sets GMT source",
    )
    collection = require_nonempty_string(
        gmt["collection"], field="gene_sets.gmt.collection"
    )
    gmt_version = require_nonempty_string(
        gmt["version"], field="gene_sets.gmt.version"
    )
    identifier_type = require_nonempty_string(
        gmt["identifier_type"], field="gene_sets.gmt.identifier_type"
    )
    if identifier_type != GENE_SET_IDENTIFIER_TYPE:
        raise InvalidRequestError(
            "'gene_sets.gmt.identifier_type' must be 'symbol'.",
            details={"observed": identifier_type},
        )

    annotation = value["annotation"]
    if not isinstance(annotation, Mapping):
        raise InvalidRequestError("'gene_sets.annotation' must be an object.")
    require_expected_keys(
        annotation,
        allowed=_ANNOTATION_KEYS,
        required=_ANNOTATION_KEYS,
        context="analysis gene_sets annotation source",
    )
    annotation_name = require_nonempty_string(
        annotation["name"], field="gene_sets.annotation.name"
    )
    annotation_version = require_nonempty_string(
        annotation["version"], field="gene_sets.annotation.version"
    )
    gene_id_column = require_nonempty_string(
        annotation["gene_id_column"], field="gene_sets.annotation.gene_id_column"
    )
    symbol_column = require_nonempty_string(
        annotation["symbol_column"], field="gene_sets.annotation.symbol_column"
    )
    if gene_id_column != GENE_SET_GENE_ID_COLUMN:
        raise InvalidRequestError(
            "'gene_sets.annotation.gene_id_column' must be 'gene_id'.",
            details={"observed": gene_id_column},
        )
    if symbol_column != GENE_SET_SYMBOL_COLUMN:
        raise InvalidRequestError(
            "'gene_sets.annotation.symbol_column' must be 'symbol'.",
            details={"observed": symbol_column},
        )

    minimum = value["minimum_tested_genes"]
    if isinstance(minimum, bool) or not isinstance(minimum, int) or minimum <= 0:
        raise InvalidRequestError(
            "'gene_sets.minimum_tested_genes' must be a positive integer."
        )

    gmt_source = _capture_source(
        gmt["path"],
        expected_sha256=gmt["sha256"],
        request_path=request_path,
        role="gmt",
    )
    annotation_source = _capture_source(
        annotation["path"],
        expected_sha256=annotation["sha256"],
        request_path=request_path,
        role="annotation",
    )
    return {
        "gmt": {
            **gmt_source,
            "collection": collection,
            "version": gmt_version,
            "identifier_type": identifier_type,
        },
        "annotation": {
            **annotation_source,
            "name": annotation_name,
            "version": annotation_version,
            "gene_id_column": gene_id_column,
            "symbol_column": symbol_column,
        },
        "minimum_tested_genes": minimum,
    }


__all__ = [
    "GENE_SET_GENE_ID_COLUMN",
    "GENE_SET_IDENTIFIER_TYPE",
    "GENE_SET_SYMBOL_COLUMN",
    "parse_gene_sets_request",
]
