"""
Main entry point for running virome_classifier as a module.

Usage:
    python -m virome_classifier [options]
"""

import sys
from .cli.classify import main

if __name__ == "__main__":
    sys.exit(main())
