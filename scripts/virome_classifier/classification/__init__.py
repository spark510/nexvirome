"""
Classification modules.

- EMClassifier: EM algorithm for multi-mapping read resolution
"""

from .em_classifier import EMClassifier
from .lca_classifier import LCAClassifier

__all__ = ["EMClassifier", "LCAClassifier"]
