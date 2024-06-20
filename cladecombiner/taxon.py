from typing import Any


class Taxon:
    """Representation of taxonomic units"""

    def __init__(self, name: str, is_tip: bool, data: Any = None):
        if not isinstance(is_tip, bool):
            raise TypeError()
        self.name: str = name
        self.tip: bool = is_tip
        self.data: Any = data

    def __eq__(self, other) -> bool:
        return self.name == other.name and self.tip == other.tip

    def __hash__(self) -> int:
        return hash(str(self.name) + str(self.tip))

    def __repr__(self) -> str:
        return f"Taxon({self.name}, tip={str(self.tip)})"
