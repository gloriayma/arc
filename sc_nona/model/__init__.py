"""Make the local NTv3 source importable for all submodules.

This is done once when the `model` package is first imported, so individual
submodules can simply do `from modeling_ntv3_pretrained import ...`.
"""
import sys
from pathlib import Path

_NTV3_DIR = Path(__file__).resolve().parents[2] / "ntv3" / "ntv3_base_model"
if str(_NTV3_DIR) not in sys.path:
    sys.path.insert(0, str(_NTV3_DIR))
