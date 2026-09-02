# Legacy HOMER extra

This directory preserves the pre-P0 promoter-motif helper for historical use.
It is not installed with `rnaseq-downstream`, is not part of the default
Snakemake DAG, and is not a certified toolkit capability. The implementation
uses a thresholded differential-gene list and HOMER's external
`findMotifs.pl`; it has no P0 statistical or benchmark evidence.

The retained runner expects the repository's legacy optional dependencies and
configuration. It is provided as unmaintained reference code only.
