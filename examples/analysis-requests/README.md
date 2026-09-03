# Analysis-request examples

`v1.1-pathways.request.json` demonstrates the complete version 1.1 shape for
the optional frozen local gene-set analysis. The bundled GMT and annotation
files are deliberately synthetic schema examples, not a biological resource.
Replace their paths, identity labels, and SHA-256 digests with a reviewed,
versioned local collection and matching annotation table.

The `validated_input_bundle` path is a placeholder. It must point to the
immutable directory produced by `rnaseq-downstream validate` for the samples
being analyzed.
