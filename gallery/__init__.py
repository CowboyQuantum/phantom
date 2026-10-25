"""The gallery renderer (ADR 0008, 0009): runs `docs/gallery/protocols/*.yaml`, hashes the arrays,
checks sign-offs and budgets, renders PNGs and an index. Not part of the `phantom` package; it imports
phantom and matplotlib, never the reverse."""
