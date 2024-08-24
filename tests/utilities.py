import os
import numpy as np
import shutil
import subprocess


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
        if (len(l) > 1) & (original_start_init is None):
            original_start_init = cpt-1
    if len(original_end_init) > 0:  
        original_end_init = original_end_init[0]
    else:
        original_end_init = last_line
    if original_end_init < original_start_init:
        original_end_init = last_line
    return original_start_init, original_end_init+1

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
    list_of_empty_line = ["\n"]
    for _ in np.arange(20):
        list_of_empty_line.append(" "+list_of_empty_line[-1])
    status = []
    for ncpt, nl in enumerate(ncontent):
        if nl in list_of_empty_line:
            status.append(["LINE IS EMPTY", ncpt, None])
        else:
            FOUND = False
            for ocpt, ol in enumerate(ocontent):
                if len(ol) > 1:
                    if nl in ol:
                        FOUND = True
                        oloc = ocpt
            if FOUND:
                status.append(["LINE DOES EXIST", ncpt, oloc])
            else:
                if "(...)" in nl:
                    status.append(["TRANSITION", ncpt, None])
                else:
                    status.append(["LINE DOESNT EXIST", ncpt, None])

    return np.array(status)

def write_non_method(new_content, original_file_content, folder, name, indent):
    # check if the new lines are new or not
    status_new_lines = detect_existing_lines(new_content,
                                             original_file_content)   
    # The new content is "import" function, to be added at the
    # beginning of the file
    file = open(folder+name+".py", "w")
    for line, status_line in zip(new_content, status_new_lines):
        status, pos_new, pos_old = status_line
        if "LINE DOESNT EXIST" in status: # only add non empty line that don't exists 
            file.write(indent+line)
    for line in original_file_content:
        file.write(line)
    file.close()

def write_new_method(new_content, original_file_content, folder, name, indent, status_new_lines):

    file = open(folder+name+".py", "w")
    for line in original_file_content:
        file.write(line)
    file.write("\n")
    for line, status_line in zip(new_content, status_new_lines):
        status, pos_new, pos_old = status_line
        if "LINE DOESNT EXIST" in status: # only add non empty line that don't exists 
            file.write(indent+line)
    file.close()

def append_new_line_to_method(new_content, original_file_content, folder, name,
                              indent, status_new_lines, original_start):

    # transition_id = np.where(status_new_lines[:,0] == 'TRANSITION')[0]
    # new_line_id = np.where(status_new_lines[:,0] == 'LINE DOESNT EXIST')[0]
    last_exising_line = np.max(status_new_lines[status_new_lines[:,0] == 'LINE DOES EXIST'][:,2])
    last_exising_line += 1
    #if True in np.unique(new_line_id > transition_id):    
    # the new line are after the transition
    # new content to be placed at the end of the method
    file = open(folder+name+".py", "w")
    for line in original_file_content[:original_start+last_exising_line]:
        file.write(line)
    for line, status_line in zip(new_content,
                                status_new_lines):
        status, pos_new, pos_old = status_line
        if "LINE DOESNT EXIST" in status: # only add non empty line that don't exists 
            file.write(indent+line)
    for line in original_file_content[original_start+last_exising_line:]:
        file.write(line)
    file.close()

def append_new_line_no_transition(new_content, original_file_content, folder, name, indent, status_new_lines, original_start):
    file = open(folder+name+".py", "w")
    for cpt, line in enumerate(original_file_content):
        file.write(line)
        if cpt-original_start in status_new_lines[:,2]:
            STOP = True
            for new_line, status_line in zip(new_content, status_new_lines):
                status, pos_new, pos_old = status_line
                try:
                    if pos_old == cpt-original_start:
                        STOP = False
                    elif pos_old > cpt-original_start:
                        STOP = True
                except:
                    assert pos_old is None
                if ("LINE DOESNT EXIST" in status) & (pos_old is None) & (STOP is False): # only add non empty line that don't exists 
                    file.write(indent+new_line)
    file.close()

def append_content(folder, name, new_content, type):
    if "test_" not in name:
        ISIMPORT, ISMETHOD, ISCLASS, ISPARTIAL = type
        indent = ""
        if np.sum(type) > 0: # nothing to append
            original_file_content = return_file_content(folder+name+".py")
            existing_methods = detect_methods(original_file_content)
            new_method = detect_methods(new_content)

            if len(new_method) == 0: # the added block is not a method
                assert ISIMPORT
                write_non_method(new_content, original_file_content,
                                 folder, name, indent)
            elif len(new_method) == 1: # the added block is a method

                # if the class is not specified in the doc,
                # the intend will be wrong
                if ISCLASS is False:
                    indent = "    "

                original_start, original_end = detect_method_boundaries(new_method[0],
                                                                        original_file_content)
                new_start, new_end = detect_method_boundaries(new_method[0],
                                                              new_content)
                # detect the line that are in both original and readded file
                status_new_lines = detect_existing_lines(new_content[new_start:new_end],
                                                         original_file_content[original_start:original_end])   
                
                NOTHING_TO_ADD = False
                if (len(np.unique(status_new_lines[:,0])) == 1):
                    if ("LINE DOES EXIST" in status_new_lines[:,0]) | ("LINE IS EMPTY" in status_new_lines[:,0]):
                        NOTHING_TO_ADD = True
                if (len(np.unique(status_new_lines[:,0])) == 2):
                    if ("LINE DOES EXIST" in status_new_lines[:,0]) & ("LINE IS EMPTY" in status_new_lines[:,0]):
                        NOTHING_TO_ADD = True

                if NOTHING_TO_ADD == False:
                    if new_method[0] not in existing_methods:
                        # this is a new method, to be added at the end
                        write_new_method(new_content[new_start:new_end], original_file_content,
                                        folder, name, indent, status_new_lines)
                    else:                            
                        if np.sum(status_new_lines[:,0] == 'TRANSITION') == 1:
                            # there is one transition (...) in new_content
                            append_new_line_to_method(new_content[new_start:new_end], original_file_content,
                                                    folder, name, indent, status_new_lines, original_start)
                        elif np.sum(status_new_lines[:,0] == 'TRANSITION') == 0:
                            # There if no transition, but content must be added to existing method
                            append_new_line_no_transition(new_content[new_start:new_end], original_file_content,
                                                            folder, name, indent, status_new_lines, original_start)
                        elif np.sum(status_new_lines[:,0] == 'TRANSITION') == 2:
                            append_new_line_to_method(new_content[new_start:new_end], original_file_content,
                                                    folder, name, indent, status_new_lines, original_start)
                        else:
                            print("NOT ANTICIPATED: SEVERAL METHODS")
            else:
                print("NOT ANTICIPATED: SEVERAL METHODS")

def sphinx_to_python(path_to_docs, chapter_id):
    filename = path_to_docs + "chapter"+str(chapter_id)+".rst"
    RST_EXISTS = os.path.exists(filename)
    created_files, created_tests = [], []
    if RST_EXISTS:  
        folder = detect_saving_folder(chapter_id)
        created_files = try_to_copy_file(chapter_id, created_files)
        file_content = return_file_content(filename)
        block_contents, block_names = detect_block_code(file_content)
        block_types = detect_block_types(block_contents)
        created_files, created_tests = create_file(block_contents, block_names,
                                                    created_files, created_tests, folder)
        for new_content, name, type in zip(block_contents, block_names, block_types):
            append_content(folder, name, new_content, type)
        return RST_EXISTS, created_tests, folder
    else:
        return RST_EXISTS, created_files, None

def run_the_test(folder, created_tests, chapter_id):
    mycwd = os.getcwd() # initial path
    os.chdir(folder)
    for test_file in created_tests:
        print("TEST --", "chapter"+str(chapter_id)+".rst", "--", test_file)
        subprocess.call(["python3", test_file])
    os.chdir(mycwd)