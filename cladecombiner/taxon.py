from typing import Any


class Taxon:
    """
    Representation of taxonomic units
    """

    def __init__(self, name: str, is_tip: bool, data: Any = None):
        self.name: str = name
        self.tip: bool = is_tip
        self.data: Any = data

        if not isinstance(is_tip, bool):
            raise TypeError(f"is_tip is {is_tip}, not a boolean")

    def __eq__(self, other) -> bool:
        return self.name == other.name and self.tip == other.tip

    def __repr__(self) -> str:
        return f"Taxon({repr(self.name)}, tip={repr(self.tip)})"
