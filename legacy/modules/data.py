import csv
import json
import math
import os
import re
from collections import Counter
from decimal import Decimal, InvalidOperation
from itertools import islice
from numbers import Integral

import pandas as pd
import yaml

from rnaseq_downstream.errors import (
    CountValuesInvalidError,
    GeneIdentifierError,
    InputReadError,
    InputValidationError,
    SampleSetMismatchError,
)


def normalize_gene_ids(index, strip_gene_version=True):
    values = pd.Index(index).astype(str)
    if strip_gene_version:
        return values.str.replace(r'\.\d+$', '', regex=True)
    return values


def _validate_and_normalize_gene_ids(
    values,
    *,
    path,
    column,
    strip_gene_version,
):
    raw = pd.Series(values, copy=True).reset_index(drop=True)
    if raw.empty:
        raise GeneIdentifierError(
            "Gene identifier column must contain at least one row.",
            details={
                'path': str(path),
                'column': column,
                'reason': 'no_gene_rows',
            },
        )

    missing = raw.isna() | raw.astype('string').str.strip().eq('').fillna(True)
    if missing.any():
        raise GeneIdentifierError(
            "Gene identifiers must be present and non-empty.",
            details={
                'path': str(path),
                'column': column,
                'reason': 'missing_or_empty',
                'row_numbers': [int(position) + 2 for position in raw.index[missing]],
            },
        )

    text = raw.astype(str)
    surrounding_whitespace = text.ne(text.str.strip())
    if surrounding_whitespace.any():
        raise GeneIdentifierError(
            "Gene identifiers must not contain surrounding whitespace.",
            details={
                'path': str(path),
                'column': column,
                'reason': 'surrounding_whitespace',
                'row_numbers': [
                    int(position) + 2 for position in text.index[surrounding_whitespace]
                ],
            },
        )

    normalized = pd.Series(
        normalize_gene_ids(text, strip_gene_version=strip_gene_version),
        index=text.index,
    )
    normalized_empty = normalized.str.strip().eq('')
    if normalized_empty.any():
        raise GeneIdentifierError(
            "Gene identifiers must remain non-empty after normalization.",
            details={
                'path': str(path),
                'column': column,
                'reason': 'empty_after_normalization',
                'strip_gene_version': bool(strip_gene_version),
                'source_gene_ids': sorted(text[normalized_empty].unique().tolist()),
            },
        )

    duplicated = normalized.duplicated(keep=False)
    if duplicated.any():
        collisions = {}
        for gene_id in sorted(normalized[duplicated].unique().tolist()):
            source_ids = sorted(text[normalized.eq(gene_id)].unique().tolist())
            collisions[str(gene_id)] = source_ids
        raise GeneIdentifierError(
            "Gene identifiers must remain unique after normalization.",
            details={
                'path': str(path),
                'column': column,
                'reason': 'duplicate_or_normalization_collision',
                'strip_gene_version': bool(strip_gene_version),
                'collisions': collisions,
            },
        )

    return pd.Index(normalized, name='gene_id')


def _read_table(path, *, operation, **read_csv_kwargs):
    try:
        return pd.read_csv(path, **read_csv_kwargs)
    except Exception as error:
        raise InputReadError(
            "Failed to read tabular input.",
            path=path,
            operation=operation,
            cause=error,
        ) from error


def _read_delimited_header(path, *, delimiter, operation):
    try:
        with open(path, 'r', encoding='utf-8', newline='') as handle:
            if delimiter is None:
                sample = handle.read(8192)
                handle.seek(0)
                try:
                    dialect = csv.Sniffer().sniff(sample, delimiters=',\t;|')
                except csv.Error:
                    first_line = handle.readline().rstrip('\r\n')
                    header = [first_line] if first_line else []
                    return header, '\t'
                return next(csv.reader(handle, dialect=dialect)), dialect.delimiter
            return next(csv.reader(handle, delimiter=delimiter)), delimiter
    except Exception as error:
        raise InputReadError(
            "Failed to read tabular input header.",
            path=path,
            operation=operation,
            cause=error,
        ) from error


def _validate_row_widths(path, header, *, delimiter, operation):
    """Reject malformed records before pandas can silently reinterpret columns."""

    invalid_rows = []
    invalid_count = 0
    try:
        with open(path, 'r', encoding='utf-8', newline='') as handle:
            if delimiter == r'\s+':
                records = (
                    (line_number, line.rstrip('\r\n').split())
                    for line_number, line in enumerate(handle, start=1)
                )
            else:
                reader = csv.reader(handle, delimiter=delimiter)
                records = ((reader.line_num, row) for row in reader)

            next(records, None)
            for row_number, row in records:
                if len(row) == len(header):
                    continue
                invalid_count += 1
                if len(invalid_rows) < 20:
                    invalid_rows.append({
                        'row_number': int(row_number),
                        'observed_field_count': len(row),
                    })
    except (OSError, UnicodeError, csv.Error) as error:
        raise InputReadError(
            "Failed to validate tabular input records.",
            path=path,
            operation=operation,
            cause=error,
        ) from error

    if invalid_count:
        raise InputValidationError(
            "Every input record must contain exactly the declared number of fields.",
            details={
                'path': str(path),
                'reason': 'row_field_count_mismatch',
                'expected_field_count': len(header),
                'invalid_row_count': invalid_count,
                'examples': invalid_rows,
            },
        )


def _read_metadata_header(path, *, operation):
    try:
        with open(path, 'r', encoding='utf-8', newline='') as handle:
            first_line = handle.readline().rstrip('\r\n')
            if '\t' in first_line:
                return next(csv.reader([first_line], delimiter='\t')), '\t'
            return first_line.split(), r'\s+'
    except Exception as error:
        raise InputReadError(
            "Failed to read tabular input header.",
            path=path,
            operation=operation,
            cause=error,
        ) from error


def _duplicate_header_occurrences(header):
    counts = Counter(header)
    return {
        column: counts[column] for column in sorted(counts) if counts[column] > 1
    }


def _validate_expression_header(path):
    header, delimiter = _read_delimited_header(
        path,
        delimiter='\t',
        operation='load_expression_matrix',
    )
    blank_positions = [
        index + 1 for index, column in enumerate(header) if not str(column).strip()
    ]
    if blank_positions:
        raise InputValidationError(
            "Expression matrix column names must be non-empty.",
            details={
                'path': str(path),
                'reason': 'empty_expression_column_name',
                'column_positions': blank_positions,
            },
        )
    occurrences = _duplicate_header_occurrences(header)
    if not occurrences:
        _validate_row_widths(
            path,
            header,
            delimiter=delimiter,
            operation='load_expression_matrix',
        )
        return header

    if 'gene_id' in occurrences:
        raise GeneIdentifierError(
            "Expression matrix must contain exactly one 'gene_id' column.",
            details={
                'path': str(path),
                'column': 'gene_id',
                'reason': 'duplicate_gene_id_column',
                'occurrences': occurrences['gene_id'],
            },
        )

    duplicate_samples = sorted(
        column for column in occurrences if column not in {'gene_id', 'gene_name'}
    )
    if duplicate_samples:
        raise SampleSetMismatchError(
            "Expression matrix sample columns must be unique.",
            details={
                'path': str(path),
                'reason': 'duplicate_matrix_samples',
                'duplicate_samples': duplicate_samples,
                'occurrences': {
                    sample: occurrences[sample] for sample in duplicate_samples
                },
            },
        )

    raise InputValidationError(
        "Expression matrix display columns must be unique.",
        details={
            'path': str(path),
            'reason': 'duplicate_display_columns',
            'duplicate_columns': sorted(occurrences),
            'occurrences': occurrences,
        },
    )


def _validate_annotation_header(path, *, delimiter, gene_id_col):
    header, delimiter = _read_delimited_header(
        path,
        delimiter=delimiter,
        operation='load_annotation_table',
    )
    blank_positions = [
        index + 1 for index, column in enumerate(header) if not str(column).strip()
    ]
    if blank_positions:
        raise InputValidationError(
            "Annotation table column names must be non-empty.",
            details={
                'path': str(path),
                'reason': 'empty_annotation_column_name',
                'column_positions': blank_positions,
            },
        )
    occurrences = _duplicate_header_occurrences(header)
    if gene_id_col in occurrences:
        raise GeneIdentifierError(
            "Annotation table must contain exactly one gene ID column.",
            details={
                'path': str(path),
                'column': gene_id_col,
                'reason': 'duplicate_gene_id_column',
                'occurrences': occurrences[gene_id_col],
            },
        )
    if occurrences:
        raise InputValidationError(
            "Annotation table column names must be unique.",
            details={
                'path': str(path),
                'reason': 'duplicate_annotation_columns',
                'duplicate_columns': sorted(occurrences),
                'occurrences': occurrences,
            },
        )
    _validate_row_widths(
        path,
        header,
        delimiter=delimiter,
        operation='load_annotation_table',
    )
    return delimiter


def _validate_metadata_header(path):
    header, separator = _read_metadata_header(path, operation='load_metadata')
    blank_positions = [
        index + 1 for index, column in enumerate(header) if not str(column).strip()
    ]
    if blank_positions:
        raise InputValidationError(
            "Metadata column names must be non-empty.",
            details={
                'path': str(path),
                'reason': 'empty_metadata_column_name',
                'column_positions': blank_positions,
            },
        )
    occurrences = _duplicate_header_occurrences(header)
    if occurrences:
        raise InputValidationError(
            "Metadata column names must be unique.",
            details={
                'path': str(path),
                'reason': 'duplicate_metadata_columns',
                'duplicate_columns': sorted(occurrences),
                'occurrences': occurrences,
            },
        )
    _validate_row_widths(
        path,
        header,
        delimiter=separator,
        operation='load_metadata',
    )
    return header, separator


def _build_gene_name_annotation(values, gene_ids):
    """Build a one-row-per-gene display mapping without using symbols as keys."""

    if values is None:
        gene_names = pd.Series(
            pd.array([pd.NA] * len(gene_ids), dtype='string'),
            index=gene_ids,
            name='gene_name',
        )
    else:
        gene_names = pd.Series(values, dtype='string').reset_index(drop=True)
        missing = gene_names.isna() | gene_names.str.strip().eq('').fillna(True)
        gene_names = gene_names.mask(missing, pd.NA)
        gene_names.index = gene_ids
        gene_names.name = 'gene_name'
    return gene_names.to_frame()


def merge_gene_name_annotations(*annotations):
    """Merge display mappings and reject conflicting symbols for a stable ID."""

    available = [annotation for annotation in annotations if annotation is not None]
    if not available:
        return None

    combined = available[0][['gene_name']].copy()
    if combined.index.has_duplicates:
        raise GeneIdentifierError(
            "Gene-name annotation indices must be unique stable gene IDs.",
            details={'reason': 'ambiguous_gene_name_mapping'},
        )

    for annotation in available[1:]:
        current = annotation[['gene_name']].copy()
        if current.index.has_duplicates:
            raise GeneIdentifierError(
                "Gene-name annotation indices must be unique stable gene IDs.",
                details={'reason': 'ambiguous_gene_name_mapping'},
            )

        overlap = combined.index.intersection(current.index, sort=False)
        left = combined.loc[overlap, 'gene_name']
        right = current.loc[overlap, 'gene_name']
        conflicts = left.notna() & right.notna() & left.astype(str).ne(right.astype(str))
        if conflicts.any():
            conflicting = overlap[conflicts.to_numpy()]
            conflict_ids = conflicting[:20].astype(str).tolist()
            examples = [
                {
                    'gene_id': str(gene_id),
                    'gene_names': sorted(
                        {
                            str(combined.at[gene_id, 'gene_name']),
                            str(current.at[gene_id, 'gene_name']),
                        }
                    ),
                }
                for gene_id in conflicting[:20]
            ]
            raise GeneIdentifierError(
                "A stable gene ID maps to conflicting display symbols.",
                details={
                    'reason': 'ambiguous_gene_name_mapping',
                    'conflicting_gene_id_count': len(conflicting),
                    'conflicting_gene_ids': conflict_ids,
                    'conflicting_gene_ids_truncated': len(conflicting) > len(conflict_ids),
                    'examples': examples,
                },
            )

        union = combined.index.union(current.index, sort=False)
        combined = combined.reindex(union)
        current = current.reindex(union)
        combined['gene_name'] = combined['gene_name'].combine_first(
            current['gene_name']
        )

    combined.index.name = 'gene_id'
    return combined


def load_annotation_table(
    path,
    gene_id_col='gene_id',
    gene_name_col='gene_name',
    strip_gene_version=True,
):
    if not path:
        return None

    suffix = os.path.splitext(os.fspath(path))[1].lower()
    if suffix in {'.tsv', '.txt'}:
        separator = '\t'
    elif suffix == '.csv':
        separator = ','
    else:
        separator = None
    separator = _validate_annotation_header(
        path,
        delimiter=separator,
        gene_id_col=gene_id_col,
    )
    df = _read_table(
        path,
        operation='load_annotation_table',
        sep=separator,
        engine='python',
        dtype={gene_id_col: 'string', gene_name_col: 'string'},
        keep_default_na=False,
    )
    if gene_id_col not in df.columns:
        raise GeneIdentifierError(
            f"Annotation file is missing gene ID column '{gene_id_col}'.",
            details={
                'path': str(path),
                'column': gene_id_col,
                'reason': 'column_missing',
                'available_columns': [str(column) for column in df.columns],
            },
        )

    gene_ids = _validate_and_normalize_gene_ids(
        df[gene_id_col],
        path=path,
        column=gene_id_col,
        strip_gene_version=strip_gene_version,
    )
    gene_names = df[gene_name_col] if gene_name_col in df.columns else None
    return _build_gene_name_annotation(gene_names, gene_ids)


def load_upstream_manifest(path):
    if not path:
        return {}
    if not os.path.exists(path):
        raise FileNotFoundError(f"Upstream manifest not found: {path}")

    suffix = os.path.splitext(path)[1].lower()
    if suffix in {'.yaml', '.yml'}:
        with open(path, 'r', encoding='utf-8') as handle:
            data = yaml.safe_load(handle) or {}
    elif suffix == '.json':
        with open(path, 'r', encoding='utf-8') as handle:
            data = json.load(handle) or {}
    else:
        df = pd.read_csv(path, sep=None, engine='python')
        if df.shape[1] < 2:
            raise ValueError(
                "Tabular upstream manifest must contain at least two columns: key and value."
            )
        data = dict(zip(df.iloc[:, 0].astype(str), df.iloc[:, 1].astype(str)))

    if not isinstance(data, dict):
        raise ValueError(f"Upstream manifest must define a mapping: {path}")
    return data


def build_upstream_provenance(config_values, manifest=None):
    manifest = manifest or {}
    provenance = {
        'upstream_manifest': config_values.get('UPSTREAM_MANIFEST'),
        'upstream_pipeline_name': config_values.get('UPSTREAM_PIPELINE_NAME'),
        'upstream_pipeline_version': config_values.get('UPSTREAM_PIPELINE_VERSION'),
        'upstream_pipeline_url': config_values.get('UPSTREAM_PIPELINE_URL'),
        'reference_genome': config_values.get('REFERENCE_GENOME'),
        'annotation_release': config_values.get('ANNOTATION_RELEASE'),
        'quantification_method': config_values.get('QUANTIFICATION_METHOD'),
        'count_matrix_type': config_values.get('COUNT_MATRIX_TYPE'),
        'tpm_matrix_type': config_values.get('TPM_MATRIX_TYPE'),
    }

    for key, value in manifest.items():
        provenance.setdefault(str(key), value)
        if provenance.get(str(key)) in [None, '']:
            provenance[str(key)] = value

    return {key: value for key, value in provenance.items() if value not in [None, '']}


def load_expression_matrix(
    path,
    metadata_samples=None,
    strip_gene_version=True,
    aggregate=None,
    validate_numeric=True,
    preserve_sample_lexemes=False,
    return_gene_annotations=False,
):
    """Load a gene-by-sample matrix without aggregating or dropping columns.

    ``gene_id`` is the required stable primary key. ``gene_name``, when present,
    is display metadata and is never used as an index. Numeric expression is
    validated without missing-value imputation; ``load_counts`` opts into its
    stricter count-specific validator instead. ``aggregate`` is retained only
    for compatibility with legacy callers; duplicates always fail instead of
    being combined, regardless of its value.
    """
    header = _validate_expression_header(path)
    sample_columns = [
        column for column in header if column not in {'gene_id', 'gene_name'}
    ]
    column_types = {'gene_id': 'string'}
    if 'gene_name' in header:
        column_types['gene_name'] = 'string'
    if preserve_sample_lexemes:
        column_types.update({column: 'string' for column in sample_columns})
    df = _read_table(
        path,
        operation='load_expression_matrix',
        sep='\t',
        dtype=column_types,
        keep_default_na=False,
    )
    if 'gene_id' not in df.columns:
        raise GeneIdentifierError(
            "Expression matrix must contain an explicit 'gene_id' column.",
            details={
                'path': str(path),
                'column': 'gene_id',
                'reason': 'column_missing',
                'available_columns': [str(column) for column in df.columns],
            },
        )

    gene_ids = _validate_and_normalize_gene_ids(
        df['gene_id'],
        path=path,
        column='gene_id',
        strip_gene_version=strip_gene_version,
    )
    if not sample_columns:
        raise InputValidationError(
            "Expression matrix must contain at least one sample column.",
            details={
                'path': str(path),
                'reason': 'no_sample_columns',
            },
        )
    matrix = df.loc[:, sample_columns].copy()
    matrix.index = gene_ids
    annotations = _build_gene_name_annotation(
        df['gene_name'] if 'gene_name' in df.columns else None,
        gene_ids,
    )
    if validate_numeric:
        matrix = _validate_numeric_expression(matrix, path=path)

    if metadata_samples is None:
        if return_gene_annotations:
            return matrix, annotations
        return matrix

    requested_samples = list(metadata_samples)
    duplicated_metadata = pd.Index(requested_samples).duplicated(keep=False)
    if duplicated_metadata.any():
        duplicates = sorted(
            pd.Index(requested_samples)[duplicated_metadata].astype(str).unique().tolist()
        )
        raise InputValidationError(
            "Metadata sample identifiers must be unique.",
            details={
                'path': str(path),
                'reason': 'duplicate_metadata_samples',
                'duplicate_samples': duplicates,
            },
        )

    metadata_set = set(requested_samples)
    matrix_set = set(sample_columns)
    missing = [sample for sample in requested_samples if sample not in matrix_set]
    extra = [sample for sample in sample_columns if sample not in metadata_set]
    if missing or extra:
        raise SampleSetMismatchError(
            "Metadata and expression matrix sample sets must match exactly.",
            details={
                'path': str(path),
                'missing_from_matrix': [str(sample) for sample in missing],
                'extra_in_matrix': [str(sample) for sample in extra],
                'metadata_sample_count': len(requested_samples),
                'matrix_sample_count': len(sample_columns),
            },
        )

    matrix = matrix.loc[:, requested_samples]
    if return_gene_annotations:
        return matrix, annotations
    return matrix


def contrast_name(contrast, default_factor):
    factor = contrast['factor']
    base = f"{contrast['treatment']}_vs_{contrast['control']}"
    if factor == default_factor:
        return contrast.get('name') or base
    return contrast.get('name') or f"{factor}__{base}"


def load_contrasts(default_factor, inline_contrasts=None, contrasts_file=None):
    if contrasts_file:
        df = pd.read_csv(contrasts_file, sep=None, engine='python')
        required = {'treatment', 'control'}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(
                f"Contrasts file must contain columns: {', '.join(sorted(required))}. "
                f"Missing: {', '.join(sorted(missing))}"
            )
        records = []
        for row in df.to_dict(orient='records'):
            factor = row.get('factor') or default_factor
            records.append({
                'factor': str(factor),
                'treatment': str(row['treatment']),
                'control': str(row['control']),
                'name': str(row['name']) if row.get('name') else None,
            })
        return records

    records = []
    for item in inline_contrasts or []:
        if isinstance(item, dict):
            records.append({
                'factor': str(item.get('factor', default_factor)),
                'treatment': str(item['treatment']),
                'control': str(item['control']),
                'name': str(item['name']) if item.get('name') else None,
            })
        else:
            treatment, control = item
            records.append({
                'factor': str(default_factor),
                'treatment': str(treatment),
                'control': str(control),
                'name': None,
            })
    return records


def prepare_metadata(metadata, reference_levels, continuous_factors=None):
    df = metadata.copy()
    for factor, ref in (reference_levels or {}).items():
        if factor not in df.columns:
            raise ValueError(f"Reference level provided for unknown metadata column '{factor}'.")
        levels = [str(value) for value in df[factor].dropna().astype(str).unique()]
        if str(ref) not in levels:
            raise ValueError(
                f"Reference level '{ref}' not found in metadata column '{factor}'. "
                f"Available groups: {', '.join(sorted(levels))}"
            )
        ordered_levels = [str(ref)] + [level for level in sorted(levels) if level != str(ref)]
        df[factor] = pd.Categorical(df[factor].astype(str), categories=ordered_levels, ordered=True)

    for factor in continuous_factors or []:
        if factor not in df.columns:
            raise ValueError(f"Continuous factor '{factor}' not found in metadata.")
        df[factor] = pd.to_numeric(df[factor], errors='raise')

    return df


def load_metadata(path, design_col):
    """Loads and validates metadata."""
    header, separator = _validate_metadata_header(path)
    df = _read_table(
        path,
        operation='load_metadata',
        sep=separator,
        engine='python',
        dtype='string',
        keep_default_na=False,
    )
    if df.empty:
        raise InputValidationError(
            "Metadata table must contain at least one sample row.",
            details={'path': str(path), 'reason': 'no_sample_rows'},
        )

    sample_id_column = 'sample_id' if 'sample_id' in df.columns else header[0]
    df = df.set_index(sample_id_column)

    if design_col not in df.columns:
        raise InputValidationError(
            f"Design column '{design_col}' not found in metadata.",
            details={
                'path': str(path),
                'reason': 'design_column_missing',
                'design_column': design_col,
                'available_columns': [str(column) for column in df.columns],
            },
        )
    sample_ids = pd.Series(df.index, dtype='string')
    missing_sample_ids = sample_ids.isna() | sample_ids.str.strip().eq('').fillna(True)
    if missing_sample_ids.any():
        raise InputValidationError(
            "Metadata sample identifiers must be present and non-empty.",
            details={'path': str(path), 'reason': 'missing_sample_id'},
        )
    padded_sample_ids = sample_ids.ne(sample_ids.str.strip())
    if padded_sample_ids.any():
        raise InputValidationError(
            "Metadata sample identifiers must not contain surrounding whitespace.",
            details={
                'path': str(path),
                'reason': 'sample_id_surrounding_whitespace',
                'row_numbers': [
                    int(position) + 2
                    for position in sample_ids.index[padded_sample_ids]
                ],
            },
        )
    if df.index.has_duplicates:
        dup_ids = df.index[df.index.duplicated(keep=False)].astype(str).unique().tolist()
        raise InputValidationError(
            "Metadata sample identifiers must be unique.",
            details={
                'path': str(path),
                'reason': 'duplicate_sample_ids',
                'duplicate_samples': sorted(dup_ids),
            },
        )
    design_values = df[design_col].astype('string')
    missing_design = (
        design_values.isna()
        | design_values.str.strip().eq('').fillna(True)
    )
    if missing_design.any():
        raise InputValidationError(
            f"Missing values found in design column '{design_col}'.",
            details={
                'path': str(path),
                'reason': 'missing_design_values',
                'design_column': design_col,
                'samples': sorted(df.index[missing_design].astype(str).tolist()),
            },
        )

    return df


def _invalid_count_details(matrix, mask, *, path, reason):
    positions = zip(*mask.to_numpy().nonzero())
    examples = []
    for row_position, column_position in islice(positions, 20):
        value = matrix.iat[int(row_position), int(column_position)]
        examples.append({
            'gene_id': str(matrix.index[int(row_position)]),
            'sample_id': str(matrix.columns[int(column_position)]),
            'value': repr(value),
        })
    return {
        'path': str(path),
        'reason': reason,
        'invalid_value_count': int(mask.to_numpy().sum()),
        'examples': examples,
    }


def _validate_numeric_expression(matrix, *, path):
    blank = matrix.apply(
        lambda column: column.map(
            lambda value: isinstance(value, str) and not value.strip()
        )
    )
    missing = matrix.isna() | blank
    if missing.any(axis=None):
        raise InputValidationError(
            "Expression matrix contains missing values.",
            details=_invalid_count_details(
                matrix,
                missing,
                path=path,
                reason='missing_expression_value',
            ),
        )

    numeric = matrix.apply(pd.to_numeric, errors='coerce')
    non_numeric = numeric.isna()
    if non_numeric.any(axis=None):
        raise InputValidationError(
            "Expression matrix contains non-numeric values.",
            details=_invalid_count_details(
                matrix,
                non_numeric,
                path=path,
                reason='non_numeric_expression_value',
            ),
        )

    finite = numeric.apply(lambda column: column.map(math.isfinite))
    if not finite.all(axis=None):
        raise InputValidationError(
            "Expression matrix contains non-finite values.",
            details=_invalid_count_details(
                matrix,
                ~finite,
                path=path,
                reason='non_finite_expression_value',
            ),
        )
    return numeric


_INTEGER_LEXEME = re.compile(r'^[0-9]+$')
_INT64_MAX_LEXEME = str(2**63 - 1)
_COUNT_ERROR_PRIORITY = (
    'missing',
    'non_numeric',
    'non_finite',
    'fractional',
    'negative',
    'non_integer_lexeme',
    'integer_out_of_range',
)
_COUNT_ERROR_MESSAGES = {
    'missing': "Count matrix contains missing values.",
    'non_numeric': "Count matrix contains non-numeric values.",
    'non_finite': "Count matrix contains non-finite values.",
    'fractional': "Count matrix contains fractional values.",
    'negative': "Count matrix contains negative values.",
    'non_integer_lexeme': (
        "Count matrix values must use plain non-negative integer syntax."
    ),
    'integer_out_of_range': (
        "Count matrix contains values outside the supported integer range."
    ),
}


def _parse_integer_count(value):
    """Parse one count exactly, rejecting any lossy or inferred representation."""

    if value is pd.NA or value is None:
        return None, 'missing'
    if isinstance(value, bool):
        return None, 'non_numeric'
    if isinstance(value, Integral):
        integer = int(value)
        if integer < 0:
            return None, 'negative'
        if integer > 2**63 - 1:
            return None, 'integer_out_of_range'
        return integer, None
    if not isinstance(value, str):
        return None, 'non_integer_lexeme'

    if not value.strip():
        return None, 'missing'
    if value != value.strip():
        return None, 'non_integer_lexeme'
    if _INTEGER_LEXEME.fullmatch(value):
        normalized = value.lstrip('0') or '0'
        if (
            len(normalized) > len(_INT64_MAX_LEXEME)
            or (
                len(normalized) == len(_INT64_MAX_LEXEME)
                and normalized > _INT64_MAX_LEXEME
            )
        ):
            return None, 'integer_out_of_range'
        return int(normalized), None

    try:
        decimal = Decimal(value)
    except InvalidOperation:
        return None, 'non_numeric'
    if not decimal.is_finite():
        return None, 'non_finite'
    if decimal < 0:
        return None, 'negative'
    if decimal != decimal.to_integral_value():
        return None, 'fractional'
    return None, 'non_integer_lexeme'


def _validate_integer_counts(matrix, *, path):
    values = matrix.to_numpy(dtype=object, copy=False)
    parsed = [[0] * matrix.shape[1] for _ in range(matrix.shape[0])]
    invalid_counts = Counter()
    invalid_examples = {reason: [] for reason in _COUNT_ERROR_PRIORITY}

    for row_position, row in enumerate(values):
        for column_position, value in enumerate(row):
            integer, reason = _parse_integer_count(value)
            if reason is None:
                parsed[row_position][column_position] = integer
                continue
            invalid_counts[reason] += 1
            examples = invalid_examples[reason]
            if len(examples) < 20:
                examples.append({
                    'gene_id': str(matrix.index[row_position]),
                    'sample_id': str(matrix.columns[column_position]),
                    'value': repr(value),
                })

    for reason in _COUNT_ERROR_PRIORITY:
        if invalid_counts[reason]:
            raise CountValuesInvalidError(
                _COUNT_ERROR_MESSAGES[reason],
                details={
                    'path': str(path),
                    'reason': reason,
                    'invalid_value_count': int(invalid_counts[reason]),
                    'examples': invalid_examples[reason],
                },
            )

    return pd.DataFrame(
        parsed,
        index=matrix.index.copy(),
        columns=matrix.columns.copy(),
        dtype='int64',
    )


def load_counts(
    path,
    metadata_samples,
    min_counts=10,
    strip_gene_version=True,
    return_gene_annotations=False,
):
    """Load integer counts and retain genes with total count >= ``min_counts``.

    ``min_counts`` is an explicit inclusive threshold on the sum across all
    samples. No missing-value imputation, rounding, or duplicate aggregation is
    performed before applying that configured filter.
    """
    try:
        loaded = load_expression_matrix(
            path,
            metadata_samples,
            strip_gene_version=strip_gene_version,
            validate_numeric=False,
            preserve_sample_lexemes=True,
            return_gene_annotations=return_gene_annotations,
        )
    except InputReadError as error:
        cause = error.cause or error
        raise InputReadError(
            "Failed to load count input.",
            path=path,
            operation="load_counts",
            cause=cause,
        ) from error

    if return_gene_annotations:
        df, annotations = loaded
    else:
        df = loaded

    if isinstance(min_counts, bool) or not isinstance(min_counts, Integral):
        raise InputValidationError(
            "min_counts must be a non-negative integer.",
            details={
                'path': str(path),
                'reason': 'invalid_min_counts',
                'min_counts': repr(min_counts),
            },
        )
    if min_counts < 0:
        raise InputValidationError(
            "min_counts must be a non-negative integer.",
            details={
                'path': str(path),
                'reason': 'invalid_min_counts',
                'min_counts': int(min_counts),
            },
        )

    counts = _validate_integer_counts(df, path=path)
    totals = counts.astype(object).sum(axis=1)
    filtered = counts.loc[totals.ge(int(min_counts))]
    counts_transposed = filtered.T
    if return_gene_annotations:
        return counts_transposed, annotations.loc[filtered.index].copy()
    return counts_transposed


def validate_analysis_inputs(metadata, design_col, reference_levels, contrasts):
    """Fails fast on invalid design settings and contrasts."""
    if design_col not in metadata.columns:
        raise ValueError(f"Design factor '{design_col}' not found in metadata.")

    for factor, ref in (reference_levels or {}).items():
        groups = set(metadata[factor].astype(str).unique())
        if str(ref) not in groups:
            raise ValueError(
                f"Reference level '{ref}' not found in metadata column '{factor}'. "
                f"Available groups: {', '.join(sorted(groups))}"
            )

    invalid = []
    for contrast in contrasts:
        factor = contrast['factor']
        if factor not in metadata.columns:
            invalid.append(contrast_name(contrast, design_col))
            continue
        groups = set(metadata[factor].astype(str).unique())
        if str(contrast['treatment']) not in groups or str(contrast['control']) not in groups:
            invalid.append(contrast_name(contrast, design_col))

    if invalid:
        raise ValueError(
            f"Contrasts contain factors or groups absent from metadata: {', '.join(invalid)}"
        )
