#!/bin/bash
set -e

# pull the last version
git pull

jupyter nbconvert --to script build-documentation.ipynb
python3 build-documentation.py
rm build-documentation.py