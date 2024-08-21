#!/usr/bin/env python
# coding: utf-8

import sys, os, git
import shutil
import numpy as np
import subprocess

current_path = os.getcwd()
git_repo = git.Repo(current_path, search_parent_directories=True)
git_path = git_repo.git.rev_parse("--show-toplevel")

path_to_docs = git_path + "/docs/source/chapters/"

# make sure the documentaiton was found
assert os.path.exists(path_to_docs), """Documentation files not found"""

if os.path.exists("generated-codes/") is False:
    os.mkdir("generated-codes/")

# In chapter 1, all the files are simply created.
chapter_id = 1
filename = path_to_docs + "chapter"+str(chapter_id)+".rst"
list_files = []
if os.path.exists(filename):
    # saving folder
    folder = "generated-codes/chapter"+str(chapter_id)+"/"
    if os.path.exists(folder) is False:
        os.mkdir(folder)
    file = open(filename, "r")
    print_file = False
    for line in file: # Loop over all the lines of the file
        if ".. label::" in line: # Detect the label "start" and label "end"
            label = line.split(".. label:: ")[1] # Look for label in the line
            if label[:6] == "start_": # Detect starting label
                class_name_i = label.split("start_")[1].split("_class")[0]
                if "test" not in class_name_i:
                    list_files.append(class_name_i+".py")
                    print_file = True
                    # create file
                    myclass = open(folder+class_name_i+".py", "w")
            elif label[:4] == "end_": # Detect ending label
                class_name_f = label.split("end_")[1].split("_class")[0]
                if "test" not in class_name_f:
                    assert class_name_f == class_name_i, """Different class closed, inconsistency in rst file?"""
                    print_file = False
                    # close file
                    myclass.close()
        else:
            if print_file: # Print the content of the label into files
                if ".. code-block::" not in line: # Ignore code block line
                    if len(line) > 1: # Remove the indentation
                        myclass.write(line[4:])
                    else:
                        myclass.write(line)

for chapter_id in np.arange(2, 10):
    filename = path_to_docs + "chapter"+str(chapter_id)+".rst"
    list_classes = []
    if os.path.exists(filename):
        # saving folder
        folder = "generated-codes/chapter"+str(chapter_id)+"/"
        previous_folder = "generated-codes/chapter"+str(chapter_id-1)+"/"
        if os.path.exists(folder) is False:
            os.mkdir(folder)
        # copy all the files from the previous chapter
        for file in list_files:
            shutil.copyfile(previous_folder+"/"+file, folder+"/"+file)
        print_file = False
        for line in open(filename, "r"): # Loop over all the lines of the file
            if ".. label::" in line: # Detect the label "start" and label "end"
                label = line.split(".. label:: ")[1] # Look for label in the line
                if label[:6] == "start_": # Detect starting label
                    class_name_i = label.split("start_")[1].split("_class")[0]
                    if "test" not in class_name_i:
                        list_classes.append(class_name_i)
                        print_file = True
                        myclass = []
                elif label[:4] == "end_": # Detect ending label
                    class_name_f = label.split("end_")[1].split("_class")[0]
                    if "test" not in class_name_f:
                        assert class_name_f == class_name_i, """Different class closed, inconsistency in rst file?"""
                        print_file = False
                        # myclass.close()
                        if (len(myclass) > 0) & ("test" not in class_name_f):
                            # detect the type of code
                            ISIMPORT = False
                            ISMETHOD = False
                            ISCLASS = False
                            ISINIT = False
                            PARTIAL = False
                            new_position_class = None
                            new_position_partial = None
                            for cpt, l in enumerate(myclass):
                                if (("import" in l) & ("as" in l)) | (("from" in l) & ("import" in l)):
                                    ISIMPORT = True
                                if ("def" in l) & ("(self" in l) & ("):" in l):
                                    ISMETHOD = True
                                if ("class" in l) & (":" in l):
                                    ISCLASS = True
                                    new_position_class = cpt
                                if ("def" in l) & ("__init__" in l):
                                    ISINIT = True
                                if ("(...)" in l):
                                    PARTIAL = True
                                    new_position_partial = cpt

                            original_class = open(folder+class_name_i+".py", "r")
                            original_content = []
                            original_position_class = None
                            original_empty_lines = []
                            for cpt, l in enumerate(original_class):
                                original_content.append(l)
                                if "__init__" in l:
                                    original_position_init = cpt
                                elif "class" in l:
                                    original_position_class = cpt
                                elif "\n" == l:
                                    original_empty_lines.append(cpt)
                            original_empty_lines = np.array(original_empty_lines)
                            if ISIMPORT:
                                new_class = open(folder+class_name_i+".py", "w")
                                for l in myclass:
                                    new_class.write(l)
                                for l in original_content:
                                    new_class.write(l)
                                new_class.close()
                            elif ISMETHOD:
                                new_class = open(folder+class_name_i+".py", "w")
                                for l in original_content:
                                    new_class.write(l)
                                for l in myclass:
                                    new_class.write("    "+l)
                                new_class.close()
                            elif (ISCLASS) & (PARTIAL is False):
                                new_class = open(folder+class_name_i+".py", "w")
                                REPLACE = False
                                REPLACED = False
                                for cpt, l in enumerate(original_content):
                                    if cpt == original_position_class: # start class
                                        REPLACE = True
                                    elif (l == "\n") & REPLACE:
                                        REPLACE = False

                                    if REPLACE:
                                        if REPLACED is False:
                                            REPLACED = True
                                            for ll in myclass:
                                                new_class.write(ll)
                                    else:
                                        new_class.write(l)
                                new_class.close()
                            elif (ISINIT) & (PARTIAL):
                                end_init = original_empty_lines[original_empty_lines
                                                                > original_position_class][0]
                                new_class = open(folder+class_name_i+".py", "w")
                                for cpt, l in enumerate(original_content):
                                    new_class.write(l)
                                    if cpt == end_init+1:
                                        for ll in myclass:
                                            # Make sure the line is not already in
                                            # before writting it
                                            # ALREADYIN = False
                                            # for lll in original_content:
                                            #     if ll in lll:
                                            #         ALREADYIN = True
                                            if ("def __init__" not in ll) & ("(...)" not in ll): #  & (ALREADYIN is False):
                                                new_class.write("    "+ll)

                                new_class.close()
        
                            # remove space in empty lines
                            original_class = open(folder+class_name_i+".py", "r")
                            original_content = []
                            for cpt, l in enumerate(original_class):
                                original_content.append(l)
                            original_class.close()
                            new_class = open(folder+class_name_i+".py", "w")
                            for l in original_content:
                                if l == "    \n":
                                    l = "\n"
                                new_class.write(l)
                            new_class.close()
            else:
                if print_file: # Print the content of the label into files
                    if ".. code-block::" not in line: # Ignore code block line
                        if len(line) > 1: # Remove the indentation
                            # myclass.write(line[4:])
                            myclass.append(line[4:])
                        else:
                            # myclass.write(line)
                            myclass.append(line)

        for class_name in np.unique(list_classes):

            # Remove doublons in init functions
            original_class = open(folder+class_name+".py", "r")
            original_start_init = None
            original_end_init = []
            original_content = []
            for cpt, l in enumerate(original_class):
                original_content.append(l)
                if ("def" in l) & ("__init__" in l):
                    original_start_init = cpt
                elif (":" in l) & ("def" in l) & ("__init__" not in l):
                    original_end_init.append(cpt)        
            original_end_init = original_end_init[0]
            original_class.close()

            # detect unique lines
            _, idx = np.unique(original_content[original_start_init:original_end_init], return_index=True)
            new_bloc = []
            for i in np.sort(idx):
                new_bloc.append(original_content[original_start_init:original_end_init][i])

            REPLACED = False
            new_class = open(folder+class_name+".py", "w")
            for cpt, l in enumerate(original_content):
                if (cpt < original_start_init) | (cpt > original_end_init-2):
                    new_class.write(l)
                else:
                    if REPLACED is False:
                        REPLACED = True
                        for ll in new_bloc:
                            new_class.write(ll)

            new_class.close()

# Detect test files
for chapter_id in np.arange(1, 10):
    filename = path_to_docs + "chapter"+str(chapter_id)+".rst"
    test_files = []
    if os.path.exists(filename):
        # saving folder
        folder = "generated-codes/chapter"+str(chapter_id)+"/"
        if os.path.exists(folder) is False:
            os.mkdir(folder)
        file = open(filename, "r")
        print_file = False
        for line in file: # Loop over all the lines of the file
            if ".. label::" in line: # Detect the label "start" and label "end"
                label = line.split(".. label:: ")[1] # Look for label in the line
                if label[:6] == "start_": # Detect starting label
                    class_name_i = label.split("start_")[1].split("_class")[0]
                    if "test" in class_name_i:
                        test_files.append(class_name_i+".py")
                        print_file = True
                        # create file
                        myclass = open(folder+class_name_i+".py", "w")
                elif label[:4] == "end_": # Detect ending label
                    class_name_f = label.split("end_")[1].split("_class")[0]
                    if "test" in class_name_f:
                        assert class_name_f == class_name_i, """Different class closed, inconsistency in rst file?"""
                        print_file = False
                        # close file
                        myclass.close()
            else:
                if print_file: # Print the content of the label into files
                    if ".. code-block::" not in line: # Ignore code block line
                        if len(line) > 1: # Remove the indentation
                            myclass.write(line[4:])
                        else:
                            myclass.write(line)

# Run test files
for chapter_id in np.arange(1, 10):
    # Test
    mycwd = os.getcwd()
    os.chdir(folder)
    for test_file in test_files:
        print("TEST", "chapter"+str(chapter_id)+".rst", test_file)
        subprocess.call(["python3", test_file])
    os.chdir(mycwd)