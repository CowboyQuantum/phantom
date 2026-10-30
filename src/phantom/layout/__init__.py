"""phantom.layout — Gate polygons, the material stack (2DEG depth, dielectrics) and the material class. Owns no physics; everything above reads it. Runtime: data only.

Module `geometry`: `Gate`, `Layer`, `Stack`, `MaterialClass`, `Layout` (JSON in and out) and the textbook
`two_dot_layout()`.
"""

from phantom.layout.geometry import Gate, Layer, Layout, MaterialClass, Stack, two_dot_layout

MODEL_VERSION = 2

__all__ = ["MODEL_VERSION", "Gate", "Layer", "Layout", "MaterialClass", "Stack", "two_dot_layout"]
