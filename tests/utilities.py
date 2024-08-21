import os
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
        ISINIT = False
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
            # look for __init__ in line
            if ("def" in line) & ("__init__" in line):
                ISINIT = True
        block_types.append([ISIMPORT, ISMETHOD, ISCLASS, ISPARTIAL, ISINIT])
    return block_types

def create_file(block_contents, block_names, created_files, folder, TEST=False):
    for content, name in zip(block_contents, block_names):
        if name+".py" not in created_files:
            if (TEST is False) & ("test_" not in name):
                created_files.append(name+".py")
                file = open(folder+name+".py", "w")
                for line in content:
                    file.write(line)
                file.close()
    return created_files

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