"""
Fingerprint generator utilities.

Provides a factory function to create RDKit fingerprint generators
and module-level convenience functions for fingerprint computation.
"""

from rdkit.Chem import rdFingerprintGenerator

_GENERATOR_FACTORIES = {
    "morgan": rdFingerprintGenerator.GetMorganGenerator,
}

_DEFAULT_PARAMS = {
    "morgan": {"radius": 3, "fpSize": 2048, "includeChirality": True},
}


def get_fingerprint_generator(
    fp_type: str = "morgan", **kwargs
) -> rdFingerprintGenerator.FingerprintGenerator64:
    """
    Create and return an RDKit fingerprint generator.

    Args:
        fp_type: type of fingerprint generator. Currently supported: 'morgan'.
        **kwargs: keyword arguments forwarded to the RDKit generator constructor.
            If not provided, default parameters for the given fp_type are used.

    Returns:
        An RDKit FingerprintGenerator64 instance.

    Raises:
        ValueError: if fp_type is not supported.
    """
    if fp_type not in _GENERATOR_FACTORIES:
        raise ValueError(
            f"Unsupported fingerprint type: {fp_type!r}. "
            f"Supported types: {list(_GENERATOR_FACTORIES.keys())}"
        )
    params = {**_DEFAULT_PARAMS[fp_type], **kwargs}
    return _GENERATOR_FACTORIES[fp_type](**params)
