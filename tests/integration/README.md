# Integration tests

These tests invoke the installed CLI in a fresh subprocess with `PYTHONPATH`
removed. They verify that every outcome is represented by a single JSON
document on stdout and a stable process exit code. Checkpoint A additionally
proves read-only inspection, successful input-only validation, warning and
artifact propagation, atomic non-overwriting bundle publication, and absence
of output after scientific validation failure. It also covers raw operating-
system argument bytes at the UTF-8 JSON boundary and explicit recovery metadata
for a bundle published before parent-directory durability could be confirmed.
