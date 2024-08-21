#!/bin/bash
set -e

# pull the last version
git pull

# build the code from documentation, and run the tests
python3 build-documentation.py
