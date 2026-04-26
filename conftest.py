"""Pytest configuration: add hbb_pipeline to sys.path."""
import sys
from pathlib import Path

# Make hbb_pipeline importable without installing the package
sys.path.insert(0, str(Path(__file__).parent))
