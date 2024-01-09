# This simple script loads all Fortran77 BLAS subroutines from the reference LAPACK implementation,
# and creates a Fortran module out of them

from enum import Enum

# Read all source files from the source folder, process them, refactor them, and put all
# subroutines/function into a module
def create_fortran_module(package_name,source_folder,out_folder):

    from datetime import date
    from platform import os

    # Parameters
    today  = date.today()
    INDENT = "     "


    # Get names
    module_name = package_name + "_blas"
    module_file = module_name + ".f90"
    module_path = os.path.join(out_folder,module_file)

    # Get list of source files
    source_files = []
    for file in os.listdir(source_folder):
        if file.endswith(".f90") or file.endswith(".f") or file.endswith(".F90"):
            source_files.append(file)
    source_files.sort()

    # Parse all source files
    fortran_functions = []
    for file in source_files:
        fortran_functions.append(parse_fortran_source(source_folder,file))

    # Create file
    fid = open(module_path,"w")

    # Header
    fid.write("module {}\n".format(module_name))
    fid.write(INDENT + "implicit none(type,external)\n")
    fid.write(INDENT + "private\n\n\n\n")

    # Public interface. Assume
    for function in fortran_functions:
        fid.write(INDENT + "public :: " + function.name + "\n")

    # Actual implementation
    fid.write("\n\n" + INDENT + "contains\n\n\n\n")




    # Close module
    fid.write("end module {}\n".format(module_name))

    fid.close()

# Enum for file sections
class Section(Enum):
    HEADER = 1
    DECLARATION = 2
    BODY = 3
    END = 4

# This class represents the contents of a 1-function/1-subroutine Fortran source file parsed from BLAS/LAPACK
class Fortran_Source:
    def __init__(self):
        self.ext = ".f90"
        self.is_free_form  = True
        self.is_function   = False
        self.is_subroutine = True
        self.name          = "NONAME"
        self.body          = ""

# Read and preprocess a Fortran line for parsing: remove comments, adjust left, and if this is a continuation
# line, read all continuation lines into it
def line_read_and_preprocess(line,is_free_form):

    import re

    processed = str(line.rstrip())

    # Remove comments
    if is_free_form:
       is_comment_line = bool(re.match(r'^\s*!', processed))
       is_continuation = bool(re.match(r'^\s*&', processed))

       # If this is a continuation line, remove all that's before the continuation character
       if is_continuation:
           processed = re.sub(r'^\s*&','',processed).strip()


    else:
       is_comment_line = bool(re.match(r'^\S\S*.*', processed))
       is_continuation = bool(re.match(r'^     \S', processed))
       # Remove continuation character
       if is_continuation:
           processed = re.sub(r'^     \S', '', processed).strip()

    return processed,is_continuation,is_comment_line


def parse_fortran_source(source_folder,file_name):

    from platform import os
    import re

    # Init empty source
    Source  = Fortran_Source()
    whereAt = Section.HEADER

    if file_name.endswith(".f") or file_name.endswith(".F") or file_name.endswith(".for") or file_name.endswith(".f77"):
       Source.is_free_form = False

    # Load whole file; split by lines
    with open(os.path.join(source_folder,file_name), 'r') as file:
        # Create an empty list to store the lines
        lines = []

        # Iterate over the lines of the file
        for line in file:
            # Remove the newline character at the end of the line
            line,is_continuation,is_comment = line_read_and_preprocess(line,Source.is_free_form)

            # Append the line to the list
            if is_continuation:
               lines[-1] = lines[-1] + line
            elif is_comment:
               # Just append this line, but ensure F90+ style comment
               line = re.sub(r'^\S', '!', line)

               lines.append(line)
            else:

               # Check what section we're in
               match whereAt:
                   case Section.HEADER:
                       # Check if a declaration is starting
                       sub_found = bool('subroutine' in line.strip().lower())
                       fun_found = bool('function' in line.strip().lower())

                       if sub_found:
                           Source.is_function = False
                           Source.is_subroutine = True

                           # Find subroutine name
                           name = bool(re.match(r'^.*subroutine\s*\S+\s*\(.*\).*$',line.strip().lower()))
                           if name:
                               strip_left = re.sub(r'^.*subroutine\s*','',line.strip().lower())
                               strip_right = re.sub(r'\(.+','',strip_left.lstrip())
                               Source.name = strip_right.rstrip()

                           whereAt = Section.DECLARATION
                       elif fun_found:
                           Source.is_function = True
                           Source.is_subroutine = False

                           # Find function name
                           name = bool(re.match(r'^.*function\s*\S+\s*\(.*\).*$',line.strip().lower()))
                           if name:
                               strip_left = re.sub(r'^.*function\s*','',line.strip().lower())
                               strip_right = re.sub(r'\(.+','',strip_left.lstrip())
                               Source.name = strip_right.rstrip()

                           whereAt = Section.DECLARATION
                   #case Section.DECLARATION:

                   #case Section.BODY:

                   #case Section.END:

               # Append this line
               lines.append(line)

    #for i in range(len(lines)):
       #print(lines[i])

    return Source



# Run script
create_fortran_module("stdlib_linalg","../assets/reference_lapack/BLAS/SRC","../src")


