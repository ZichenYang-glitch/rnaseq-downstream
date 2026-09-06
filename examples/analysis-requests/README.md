# Analysis-request examples

`v1.1-pathways.request.json` demonstrates the complete version 1.1 shape for
the optional frozen local gene-set analysis. The bundled GMT and annotation
files are deliberately synthetic schema examples, not a biological resource.
Replace their paths, identity labels, and SHA-256 digests with a reviewed,
versioned local collection and matching annotation table.

`v1.2-deseq2-wald.request.json` demonstrates a coefficient-compatible apeglm
Wald request. `v1.2-deseq2-lrt.request.json` demonstrates the one reporting
contrast plus proper reduced-design subset required for an omnibus LRT.
`v1.2-edger-display.request.json` shows that version 1.2 can explicitly select
the existing edgeR display path. Omitting `backend` from that example would
still select `edger_ql`; it is included to make the routing choice visible.

The `validated_input_bundle` path is a placeholder. It must point to the
immutable directory produced by `rnaseq-downstream validate` for the samples
being analyzed.
