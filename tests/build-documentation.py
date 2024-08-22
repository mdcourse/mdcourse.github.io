
import sys, os, git
import shutil
import numpy as np
import subprocess

from utilities import detect_saving_folder, try_to_copy_file, \
                      return_file_content, detect_block_code, \
                      detect_block_types, create_file, detect_last_matching_line, \
                      append_content

current_path = os.getcwd()
git_repo = git.Repo(current_path, search_parent_directories=True)
git_path = git_repo.git.rev_parse("--show-toplevel")

path_to_docs = git_path + "/docs/source/chapters/"

# make sure the documentaiton was found
assert os.path.exists(path_to_docs), """Documentation files not found"""

if os.path.exists("generated-codes/") is False:
    os.mkdir("generated-codes/")

mycwd = os.getcwd()
for chapter_id in np.arange(10):
    filename = path_to_docs + "chapter"+str(chapter_id)+".rst"
    created_files, created_tests = [], []
    if os.path.exists(filename):
        folder = detect_saving_folder(chapter_id)
        created_files = try_to_copy_file(chapter_id, created_files)
        file_content = return_file_content(filename)
        block_contents, block_names = detect_block_code(file_content)
        block_types = detect_block_types(block_contents)
        created_files, created_tests = create_file(block_contents, block_names,
                                                   created_files, created_tests, folder)
        for content, name, type in zip(block_contents, block_names, block_types):
            append_content(folder, name, content, type)
        # Run the tests
        os.chdir(folder)
        for test_file in created_tests:
            print("TEST --", "chapter"+str(chapter_id)+".rst", "--", test_file)
            subprocess.call(["python3", test_file])
        os.chdir(mycwd)