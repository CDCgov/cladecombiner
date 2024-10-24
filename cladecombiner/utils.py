from collections.abc import Iterable
from typing import Any

def validate_kwargs(exp_kwargs: Iterable[tuple[str, Any, bool]], **kwargs):
    """
    Check that required kwargs are present and present kwargs are correct type
    """
    for arg in exp_kwargs:
        if arg[2] and (arg[0] not in kwargs):
            raise RuntimeError(f"Requires argument `{arg[0]}` is missing.")
    for arg in exp_kwargs:
        if arg[0] not in kwargs and not isinstance(kwargs[arg[0]], arg[1]):
            raise TypeError(f"Argument {arg[0]} must be of type {arg[1]}")
