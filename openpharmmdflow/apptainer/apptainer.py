"""
Module openpharmmdflow/apptainer/apptainer.py
"""

import pooch

from openpharmmdflow._version import version

GOODBOY = pooch.create(
    path=pooch.os_cache("openpharmmdflow"),
    base_url="https://github.com/omsf/OpenPharmMDFlow/releases/download/{version}/",
    version=version,
    version_dev="v0.0.1",  # Change this to 'main' once we push the tag
    env="OPENPHARMMDFLOW_DATA_DIR",
    registry={
        "cod-tools.sif": "a845ce9542931e16200a08ad172c04a2f087fcb626cdbfe3833c7993f06028cf"
    },
)
