"""
Logging system for virome_classifier.

Provides a clean, structured logging interface with different verbosity levels.
"""

import logging
import sys
from typing import Optional
from pathlib import Path
import datetime


class ViromeLogger:
    """
    Centralized logger for virome_classifier.

    Features:
    - Multiple log levels (INFO, DEBUG, VERBOSE)
    - Console and file output
    - Structured formatting
    - Performance timing
    """

    def __init__(self, name: str = "virome_classifier", verbose: bool = False):
        """
        Initialize logger.

        Args:
            name: Logger name
            verbose: Enable verbose (DEBUG) mode
        """
        self.name = name
        self._verbose = verbose
        self._logger = logging.getLogger(name)
        self._setup_logger()
        self._log_file: Optional[Path] = None

    def _setup_logger(self):
        """Setup logging configuration."""
        # Clear existing handlers
        self._logger.handlers.clear()

        # Set level
        level = logging.DEBUG if self._verbose else logging.INFO
        self._logger.setLevel(level)

        # Console handler with clean format
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(level)

        # Format: [INFO] Message
        formatter = logging.Formatter('[%(levelname)s] %(message)s')
        console_handler.setFormatter(formatter)

        self._logger.addHandler(console_handler)

        # Prevent propagation to root logger
        self._logger.propagate = False

    def set_verbose(self, verbose: bool):
        """Enable/disable verbose mode."""
        self._verbose = verbose
        self._setup_logger()

    def add_file_handler(self, log_file: Path):
        """
        Add file output handler.

        Args:
            log_file: Path to log file
        """
        self._log_file = log_file

        # Create file handler
        file_handler = logging.FileHandler(log_file, mode='a', encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)  # Always DEBUG for file

        # Detailed format for file
        formatter = logging.Formatter(
            '%(asctime)s - [%(levelname)s] %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        file_handler.setFormatter(formatter)

        self._logger.addHandler(file_handler)
        self.info(f"📝 Logging to file: {log_file}")

    def info(self, message: str):
        """Log info message."""
        self._logger.info(message)

    def debug(self, message: str):
        """Log debug message (only if verbose)."""
        self._logger.debug(message)

    def verbose(self, message: str):
        """Log verbose message (alias for debug)."""
        if self._verbose:
            self._logger.debug(f"[VERBOSE] {message}")

    def warning(self, message: str):
        """Log warning message."""
        self._logger.warning(message)

    def error(self, message: str):
        """Log error message."""
        self._logger.error(message)

    def success(self, message: str):
        """Log success message (info with ✅ prefix)."""
        self._logger.info(f"✅ {message}")

    def step(self, step_num: int, total: int, description: str):
        """
        Log pipeline step.

        Args:
            step_num: Current step number
            total: Total number of steps
            description: Step description
        """
        self._logger.info(f"{'='*80}")
        self._logger.info(f"Step {step_num}/{total}: {description}")
        self._logger.info(f"{'='*80}")


# Global logger instance
_global_logger: Optional[ViromeLogger] = None


def get_logger(name: str = "virome_classifier", verbose: bool = False) -> ViromeLogger:
    """
    Get global logger instance.

    Args:
        name: Logger name
        verbose: Enable verbose mode

    Returns:
        ViromeLogger instance
    """
    global _global_logger

    if _global_logger is None:
        _global_logger = ViromeLogger(name, verbose)

    return _global_logger


def set_verbose(verbose: bool):
    """Set verbose mode for global logger."""
    logger = get_logger()
    logger.set_verbose(verbose)


# Convenience functions for backward compatibility
def log_info(message: str):
    """Log info message."""
    get_logger().info(message)


def log_verbose(message: str):
    """Log verbose message."""
    get_logger().verbose(message)


def log_debug(message: str):
    """Log debug message."""
    get_logger().debug(message)


def log_warning(message: str):
    """Log warning message."""
    get_logger().warning(message)


def log_error(message: str):
    """Log error message."""
    get_logger().error(message)


def log_success(message: str):
    """Log success message."""
    get_logger().success(message)
