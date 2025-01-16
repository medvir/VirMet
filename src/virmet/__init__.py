"""VirMet package."""

from importlib import metadata
from typing import Final

APPLICATION_NAME: Final[str] = "virmet"
VERSION: Final[str] = metadata.version(APPLICATION_NAME)

__version__: Final[str] = VERSION
