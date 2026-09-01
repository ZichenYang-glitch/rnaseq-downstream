"""Public package surface for the RNA-seq downstream toolkit."""

from .errors import ErrorCode, ExitCode, ToolkitError

__all__ = ["ErrorCode", "ExitCode", "ToolkitError", "__version__"]

__version__ = "0.1.0.dev0"
