"""
Stub function and module used as a setuptools entry point.
"""

import augur
from sys import argv, exit
import logging

logging.basicConfig(
    format="%(levelname)-8s | %(asctime)s | %(message)s",
    level=logging.DEBUG
)


# Entry point for setuptools-installed script and bin/augur dev wrapper.
def main():
    return augur.run( argv[1:] )

# Run when called as `python -m augur`, here for good measure.
if __name__ == "__main__":
    exit( main() )
