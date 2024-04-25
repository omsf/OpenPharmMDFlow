# SPDX-FileCopyrightText: 2024-present Mike Henry <11765982+mikemhenry@users.noreply.github.com>
#
# SPDX-License-Identifier: MIT
from .io.cif import from_cif

__all__ = ["from_cif"]

from . import _version

# Add a "v" to the version number
__version__ = f"v{_version.version}"
