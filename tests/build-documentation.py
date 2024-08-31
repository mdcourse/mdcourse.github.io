import os, git
import numpy as np

# import the python converter
from utilities import sphinx_to_python, run_the_test

# detect the path to the documentation
current_path = os.getcwd()
git_repo = git.Repo(current_path, search_parent_directories=True)
git_path = git_repo.git.rev_parse("--show-toplevel")
path_to_docs = git_path + "/docs/source/chapters/"

# make sure the documentation was found
assert os.path.exists(path_to_docs), """Documentation files not found"""

# if necessary, create the "generated-codes/" folder
if os.path.exists("generated-codes/") is False:
    os.mkdir("generated-codes/")

# loop on the different chapter
for chapter_id in [1, 2, 3, 4, 5]:
    # for each chapter, create the corresponding code
    RST_EXISTS, created_tests, folder = sphinx_to_python(path_to_docs, chapter_id)
    if RST_EXISTS:
        # run the tests at the bottom of the chapter pages
        run_the_test(folder, created_tests, chapter_id)
