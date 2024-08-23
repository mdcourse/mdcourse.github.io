import os
import numpy as np
import shutil

def detect_saving_folder(chapter_id, CREATE=True):
    folder = "generated-codes/chapter"+str(chapter_id)+"/"
    if CREATE:
        if os.path.exists(folder) is False:
            os.mkdir(folder)
    return folder

def return_file_content(filename):
    file = open(filename, "r")
    file_content = []
    for line in file:
        file_content.append(line)
    file.close()
    return file_content

def detect_block_code(file_content):
    block_contents = []
    block_names = []
    LABEL = False
    for line in file_content:
        if ".. label::" in line: # Detect the label "start" and label "end"
            label = line.split(".. label:: ")[1] # Look for label in the line
            if label[:6] == "start_": # Detect starting label
                class_name_i = label.split("start_")[1].split("_class")[0]
                block = []
                LABEL = True
            elif label[:4] == "end_": # Detect ending label
                class_name_f = label.split("end_")[1].split("_class")[0]
                LABEL = False
                block_contents.append(block)
                assert class_name_i == class_name_f
                block_names.append(class_name_i)
        else:
            if LABEL: # Print the content of the label into files
                if ".. code-block::" not in line: # Ignore code block line
                    if len(line) > 1: # Remove the indentation
                        block.append(line[4:])
                    else:
                        block.append(line)
    return block_contents, block_names

def detect_block_types(block_contents):
    block_types = []
    for block in block_contents:
        ISIMPORT = False
        ISMETHOD = False
        ISCLASS = False
        ISPARTIAL = False
        for line in block:
            # look for "from X import Y" command
            if (("import" in line) & ("as" in line)) | (("from" in line) & ("import" in line)):
                ISIMPORT = True
            # look for "def X():" line
            if ("def" in line) & ("(self" in line) & ("):" in line):
                ISMETHOD = True
            # look for "class C():" line
            if ("class" in line) & (":" in line):
                ISCLASS = True
            # look for "(...)" line
            if ("(...)" in line):
                ISPARTIAL = True
        block_types.append([ISIMPORT, ISMETHOD, ISCLASS, ISPARTIAL])
    return block_types

def create_file(block_contents, block_names, created_files, created_tests, folder):
    for content, name in zip(block_contents, block_names):
        if name+".py" not in created_files:
            if "test_" not in name:
                created_files.append(name+".py")
                file = open(folder+name+".py", "w")
                for line in content:
                    file.write(line)
                file.close()
            else:
                created_tests.append(name+".py")
                file = open(folder+name+".py", "w")
                for line in content:
                    file.write(line)
                file.close()            
    return created_files, created_tests

def try_to_copy_file(chapter_id, created_files):
    new_folder = detect_saving_folder(chapter_id, CREATE=False)
    old_folder = detect_saving_folder(chapter_id-1, CREATE=False)
    if os.path.exists(old_folder):
        existing_files = next(os.walk(old_folder), (None, None, []))[2]
        for file in existing_files:
            if (".py" in file) and ("test_" not in file):
                shutil.copyfile(old_folder+"/"+file,
                                new_folder+"/"+file)
                created_files.append(file)
    return created_files

def detect_last_matching_line(content, original_file_content):
    for cpt_new, new_line in enumerate(content):
        if len(new_line) > 1:
            for cpt_old, older_line in enumerate(original_file_content):
                if len(older_line) > 0:
                    if new_line in older_line:
                        cpt_new_last = cpt_new
                        cpt_old_last = cpt_old
    return cpt_new_last, cpt_old_last

def detect_method_boundaries(method_name, original_file_content):
    original_start_init = None
    original_end_init = []
    last_line = None
    for cpt, l in enumerate(original_file_content):
        if ("def" in l) & (method_name in l):
            original_start_init = cpt
        elif (":" in l) & ("def" in l) & (method_name not in l):
            original_end_init.append(cpt)     
        if len(l) > 1:
            last_line = cpt
    if len(original_end_init) > 0:  
        original_end_init = original_end_init[0]
    else:
        original_end_init = last_line
    if original_end_init < original_start_init:
        original_end_init = last_line
    return original_start_init, original_end_init

def detect_unique_lines(file_content, start, end):
    _, idx = np.unique(file_content[start:end],
                       return_index=True)
    unique_lines = []
    for i in np.sort(idx):
        unique_lines.append(file_content[start:end][i])
    return unique_lines

def replace_method(folder, name, original_content, unique_lines, original_start, original_end):
    REPLACED = False
    new_class = open(folder+name+".py", "w")
    for cpt, l in enumerate(original_content):
        if (cpt < original_start) | (cpt > original_end-2):
            new_class.write(l)
        else:
            if REPLACED is False:
                REPLACED = True
                for ll in unique_lines:
                    new_class.write(ll)
    new_class.close()

def detect_methods(file_content):
    existing_methods = []
    for line in file_content:
        if ("def " in line) & ("(self" in line):
            method_name = line.split("def ")[1].split("(self")[0]
            existing_methods.append(method_name)
    return existing_methods

def detect_existing_lines(ncontent, ocontent):
    locations = []
    for nl in ncontent:
        if len(nl) > 1:
            FOUND = False
            for ocpt, ol in enumerate(ocontent):
                if len(ol) > 1:
                    if nl in ol:
                        FOUND = True
                        locations.append(ocpt)
        if (FOUND is False) | (len(nl) <= 1):
            # the line is new
            locations.append(-1) 
    if len(ncontent) != len(locations):
        print("WARNING")
        print(len(ncontent), len(locations))
    return locations

def append_content(folder, name, content, type):
    if "test_" not in name:
        ISIMPORT, ISMETHOD, ISCLASS, ISPARTIAL = type
        if np.sum(type) == 0:
            # nothing to append
            pass
        else:
            original_file_content = return_file_content(folder+name+".py")
            if (ISIMPORT) & (ISCLASS is False) & (ISMETHOD is False):
                # Add the content at the start of the file
                file = open(folder+name+".py", "w")
                for line in content:
                    file.write(line)
                for line in original_file_content:
                    file.write(line)
                file.close()
            else:
                existing_methods = detect_methods(original_file_content)
                new_method = detect_methods(content)
                if new_method[0] not in existing_methods:
                    # the method does not exists, it will be appened
                    assert ISPARTIAL is False
                    assert ISMETHOD
                    # Add the content at the end of the file, with an indentation
                    file = open(folder+name+".py", "w")
                    for line in original_file_content:
                        file.write(line)
                    for line in content:
                        file.write("    "+line)
                    file.close()
                elif new_method[0] in existing_methods:
                    assert (ISMETHOD) | (ISCLASS) | (ISPARTIAL)
                    original_start, original_end = detect_method_boundaries(new_method[0], original_file_content)
                    new_start, new_end = detect_method_boundaries(new_method[0], content)
                    # detect the line that are in both original and readded file
                    location_lines= detect_existing_lines(content[new_start:new_end+1],
                                                        original_file_content[original_start:original_end])

                    # if the class is not specified in the doc,
                    # the intend will be wrong
                    if ISCLASS is False:
                        indent = "    "
                    else:
                        indent = ""

                    # Add the missing lines
                    file = open(folder+name+".py", "w")
                    to_be_added = []
                    # loop on the original file
                    for ocpt, line in enumerate(original_file_content):
                        file.write(line) # add all the lines from the original file
                        if (ocpt-original_start in location_lines) & (ocpt-original_start > -1):
                            # this line exists in both the original and the new content
                            to_be_added = []
                            TOADD = False
                            ncpt = 0
                            for line, location in zip(content[new_start:new_end+1], location_lines):
                                if location == ocpt-original_start:
                                    TOADD = True
                                elif location > ocpt-original_start:
                                    TOADD = False
                                if TOADD:
                                    #if len(line) > 1:
                                    if location == -1:
                                        to_be_added.append(line)
                                ncpt += 1
                        if (ocpt-original_end in location_lines) & (ocpt-original_end > -1) & (len(to_be_added) > 0):
                            for new_line in to_be_added:
                                if "(...)" not in new_line:
                                    file.write(indent+new_line)
                    file.close()

                    # make sure there is no doublons
                    # original_file_content = return_file_content(folder+name+".py")
                    # original_start, original_end = detect_method_boundaries(new_method[0], original_file_content)
                    # unique_line = detect_unique_lines(original_file_content, original_start, original_end)
                    # file = open(folder+name+".py", "w")
                    # for cpt, line in enumerate(original_file_content):
                    #     if cpt <= original_start:
                    #         file.write(line)
                    # for line in unique_line:
                    #     file.write(line)
                    # for cpt, line in enumerate(original_file_content):
                    #     if cpt >= original_end:
                    #         file.write(line)

