# Integration tests

These tests invoke the installed CLI in a fresh subprocess with `PYTHONPATH`
removed. They verify that every outcome is represented by a single JSON
document on stdout and a stable process exit code.
