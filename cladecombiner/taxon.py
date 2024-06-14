from typing import Any


class Taxon:
    """Representation of objects to be aggregated"""

    def __init__(self, name: str, data: Any = None):
        self.name: str = name
        self.data: Any = data

    def __repr__(self) -> str:
        return f"Taxon({self.name})"
