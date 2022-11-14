"""Model the diffusion of heat over a 2D plate."""
from ._version import __version__
from .bmi_cosmic_erosion import BmiCosmicErosion
from .cosmic_erosion import CosmicErosion

__all__ = ["__version__", "BmiPaleoerosion", "Paleoerosion"]
