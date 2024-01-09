# This simple script loads all Fortran77 BLAS subroutines from the reference LAPACK implementation,
# and creates a Fortran module out of them

from enum import Enum

# Read all source files from the source folder, process them, refactor them, and put all
# subroutines/function into a module
def create_fortran_module(module_name,source_folder,out_folder,prefix):

    from datetime import date
    from platform import os

    # Parameters
    today  = date.today()
    INDENT = "     "

    # Get names
    module_file = module_name + ".f90"
    module_path = os.path.join(out_folder,module_file)

    # Get list of source files
    source_files = []
    for file in os.listdir(source_folder):
        if file.endswith(".f"): #if file.endswith(".f90") or file.endswith(".f") or file.endswith(".F90"):
            source_files.append(file)
    source_files.sort()

    # Parse all source files
    fortran_functions = []
    for file in source_files:
        fortran_functions.append(parse_fortran_source(source_folder,file,prefix))

    # Rename all procedures
    for function in fortran_functions:
        function.body,function.deps = rename_source_body(function.body,fortran_functions)

    # Create file
    fid = open(module_path,"w")

    # Header
    fid.write("module {}\n".format(module_name))
    #fid.write(INDENT + "use stdlib_kinds, only: sp,dp,lk,int32,int64\n")
    fid.write(INDENT + "use iso_fortran_env, only: int32,int64\n")
    fid.write(INDENT + "implicit none(type,external)\n")
    fid.write(INDENT + "private\n\n\n\n")

    # Temporary: to be replaced with stdlib_kinds
    fid.write(INDENT + "integer, parameter :: sp = selected_real_kind(6)\n")
    fid.write(INDENT + "integer, parameter :: dp = selected_real_kind(15)\n")
    fid.write(INDENT + "integer, parameter :: lk = kind(.true.)\n")


    # Public interface.
    fid.write(INDENT + "public :: sp,dp,lk,int32,int64\n")
    for function in fortran_functions:
        fid.write(INDENT + "public :: " + function.new_name + "\n")

        if function.new_name=="NONAME":
            print("\n".join(function.body[1:]))
            exit(1)


    # Actual implementation
    fid.write("\n\n" + INDENT + "contains\n")

    # Write functions
    print_function_tree(fortran_functions,fid)

    # Close module
    fid.write("\n\n\nend module {}\n".format(module_name))

    fid.close()

# Enum for file sections
class Section(Enum):
    HEADER = 1
    DECLARATION = 2
    EXTERNALS = 3
    BODY = 3
    END = 4

# Print function tree in a dependency-suitable way
def print_function_tree(functions,fid):

    # list of function names
    fun_names = []

    # Cleanup first
    for i in range(len(functions)):
        fun_names.append(functions[i].old_name.lower())
        functions[i].printed = False
        functions[i].ideps = []

    # Get dependency indices
    for i in range(len(functions)):
        for j in range(len(functions[i].deps)):
           functions[i].ideps.append(fun_names.index(functions[i].deps[j].lower()))
        print(functions[i].ideps)

    attempt = 0
    while attempt<len(functions):

        for i in range(len(functions)):

            # Check deps
            if (not functions[i].printed):
                nprinted = 0
                for j in range(len(functions[i].deps)):
                    # Self function: do not consider
                    if functions[i].ideps[j]==i:
                        nprinted+=1
                    elif functions[functions[i].ideps[j]].printed:
                        nprinted+=1

                if nprinted==len(functions[i].deps):
                   print(str(nprinted) + "deps printed already for " + functions[i].old_name)
                   fid.write("\n".join(functions[i].body[1:]))
                   functions[i].printed = True

        attempt+=1

    # Final check
    not_printed = 0
    for i in range(len(functions)):
        if not functions[i].printed: not_printed+=1

    if not_printed>0:
        print("***ERROR*** there are non printed functions")
        exit(1)


# This class represents the contents of a 1-function/1-subroutine Fortran source file parsed from BLAS/LAPACK
class Fortran_Source:
    def __init__(self):
        self.ext = ".f90"
        self.is_free_form  = True
        self.is_function   = False
        self.is_subroutine = True
        self.old_name      = "NONAME"
        self.new_name      = "NONAME"
        self.body          = ""
        self.deps          = []
        self.ideps         = []
        self.printed       = False

# Read and preprocess a Fortran line for parsing: remove comments, adjust left, and if this is a continuation
# line, read all continuation lines into it
def line_read_and_preprocess(line,is_free_form):

    import re

    processed = replace_f77_types(line)

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

# Parse a line, and replace old-style Fortran datatypes and constructs with stdlib kinds
def replace_f77_types(line):

    new_line = line.rstrip()
    new_line = new_line.replace(".LE.","<=")
    new_line = new_line.replace(".GE.",">=")
    new_line = new_line.replace(".EQ.","==")
    new_line = new_line.replace(".NE.","/=")
    new_line = new_line.replace(".LT.","<")
    new_line = new_line.replace(".GT.",">")
    new_line = new_line.replace("COMPLEX*16","COMPLEX(dp)")
    new_line = new_line.replace("COMPLEX ","COMPLEX(sp) ")
    new_line = new_line.replace("INTEGER ","INTEGER(int32) ")
    new_line = new_line.replace("LOGICAL ","LOGICAL(lk) ")
    new_line = new_line.replace("REAL ","REAL(sp) ")
    new_line = new_line.replace("DOUBLE PRECISION ","REAL(dp) ")
    new_line = new_line.replace("D+0","_dp")
    new_line = new_line.replace("E+0","_sp")

    return new_line

# Check if a line is a datatype line
def is_datatype_line(line):

    check_line = line.strip().lower()

    is_data_line = check_line.starts

# Check if a line is followed by external declaration
def is_externals_header(line):

    import re

    check_line = line.strip()

    # Begins with a data type
    ext =    bool(re.match(r'\S*\s*.. External Functions ..',check_line)) \
          or bool(re.match(r'\S*\s*.. External Subroutines ..',check_line))

    return ext

# Check if a line is a declaration line
def is_declaration_line(line):

    check_line = line.strip().lower()

    # Begins with a data type
    is_decl =    check_line.startswith("type ") \
              or check_line.startswith("type(") \
              or check_line.startswith("real ") \
              or check_line.startswith("real(") \
              or check_line.startswith("real::") \
              or check_line.startswith("double precision ") \
              or check_line.startswith("double precision(") \
              or check_line.startswith("double precision::") \
              or check_line.startswith("doubleprecision") \
              or check_line.startswith("integer ") \
              or check_line.startswith("integer(") \
              or check_line.startswith("integer::") \
              or check_line.startswith("complex ") \
              or check_line.startswith("complex(") \
              or check_line.startswith("complex::") \
              or check_line.startswith("character ") \
              or check_line.startswith("character(") \
              or check_line.startswith("character::") \
              or check_line.startswith("logical ") \
              or check_line.startswith("logical(") \
              or check_line.startswith("logical::") \
              or check_line.startswith("use ") \
              or check_line.startswith("use,") \
              or check_line.startswith("use::") \
              or check_line.startswith("intrinsic ") \
              or check_line.startswith("intrinsic::") \
              or check_line.startswith("external ") \
              or check_line.startswith("external::") \
              or check_line.startswith("parameter ") \
              or check_line.startswith("parameter(") \
              or check_line.startswith("parameter::")

    return is_decl

def filter_declaration_line(line):

    check_line = line.strip().lower()

    # Remove all EXTERNAL declarations
    filtered =   check_line.startswith("external ") \
              or check_line.startswith("external::")

    return filtered

# Given the list of all sources, rename all matching names in the current source body
def rename_source_body(lines,Sources):

    prefixes = [" ","(","."]

    body = []

    # Assemble list of old and new names
    old_names = []
    new_names = []
    for Source in Sources:
        old_names.append(Source.old_name.lower())
        new_names.append(Source.new_name.lower())

    is_found = [False for i in range(len(new_names))]

    for i in range(len(lines)):
        old_line = lines[i].lower()
        for j in range(len(old_names)):
            line = old_line
            for k in range(len(prefixes)):
               line = line.replace(prefixes[k] + old_names[j],prefixes[k] + new_names[j])
            if len(line)>len(old_line):
                is_found[j] = True
            old_line = line

        body.append(line);

    # Build dependency list
    dependency_list = []
    for j in range(len(old_names)):
        if is_found[j]:
            dependency_list.append(old_names[j])

    return body,dependency_list

def parse_fortran_source(source_folder,file_name,prefix):

    from platform import os
    import re

    INDENT = "     "

    # Init empty source
    Source  = Fortran_Source()
    whereAt = Section.HEADER

    if file_name.endswith(".f") or file_name.endswith(".F") or file_name.endswith(".for") or file_name.endswith(".f77"):
       Source.is_free_form = False

    # FiLoad whole file; split by lines; join concatenation lines
    with open(os.path.join(source_folder,file_name), 'r') as file:
        # Create an empty list to store the lines
        file_body = []

        # Iterate over the lines of the file
        for line in file:
            # Remove the newline character at the end of the line
            line,is_continuation,is_comment = line_read_and_preprocess(line,Source.is_free_form)

            # Append the line to the list
            if is_continuation:
               file_body[-1] = file_body[-1] + line
            else:
               file_body.append(line)

        # Iterate over the joined lines of the file
        Source.body = []

        for line in file_body:
            # Remove the newline character at the end of the line
            line,is_continuation,is_comment = line_read_and_preprocess(line,Source.is_free_form)

            # Append the line to the list
            if is_comment:

               # Inside an externals section: remove altogether
               if whereAt==Section.EXTERNALS:
                  whereAt = Section.DECLARATION
               # Is an Externals section starting
               elif whereAt==Section.DECLARATION and is_externals_header(line):
                  whereAt = Section.EXTERNALS
                  line = ""
               else:
                  # Just append this line, but ensure F90+ style comment
                  line = re.sub(r'^\S', '!', line)
                  Source.body.append(INDENT + line)

               #print("comment: "+ line.strip().lower())
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
                           name = bool(re.match(r'^.*subroutine\s*\S+\s*\(.*$',line.strip().lower()))
                           if name:
                               strip_left = re.sub(r'^.*subroutine\s*','',line.strip().lower())
                               strip_right = re.sub(r'\(.+','',strip_left.lstrip())
                               Source.old_name = strip_right.strip()
                               Source.new_name = prefix + Source.old_name

                           whereAt = Section.DECLARATION
                       elif fun_found:
                           Source.is_function = True
                           Source.is_subroutine = False

                           # Find function name
                           name = bool(re.match(r'^.*function\s*\S+\s*\(.*$',line.strip().lower()))
                           if name:
                               strip_left = re.sub(r'^.*function\s*','',line.strip().lower())
                               strip_right = re.sub(r'\(.+','',strip_left.lstrip())
                               Source.old_name = strip_right.strip()
                               Source.new_name = prefix + Source.old_name

                           whereAt = Section.DECLARATION

                   case Section.DECLARATION:

                       # A procedure name must have been read
                       if Source.new_name=="NONAME":
                           exit(1)

                       # Check if this line still begins with a declaration
                       if is_declaration_line(line):
                           # Filter declaration line
                           if (filter_declaration_line(line)):
                               line = "";
                       else:
                           # Start body section
                           whereAt = Section.BODY

                   case Section.EXTERNALS:

                       # Check if this line still begins with a declaration
                       if is_declaration_line(line):
                           # Delete altoghether
                           line = "";
                       else:
                           # Go back to declaration section
                           whereAt = Section.DECLARATION

                   case Section.BODY:

                       # End of the function/subroutine: inside a module, it must contain its name
                       if line.strip().upper()=="END" or line.strip().upper()=="END SUBROUTINE" \
                          or line.strip().upper()=="END FUNCTION":
                           whereAt = Section.END
                           if Source.is_function:
                               line = "END FUNCTION " + Source.old_name.upper()
                           elif Source.is_subroutine:
                               line = "END SUBROUTINE " + Source.old_name.upper()

                   #case Section.END:

               # Append this line
               Source.body.append(INDENT + line)

#    for i in range(len(Source.body)):
#       print(Source.body[i])

    return Source



# Run script
create_fortran_module("stdlib_linalg_blas","../assets/reference_lapack/BLAS/SRC","../src","stdlib_")
#create_fortran_module("stdlib_linalg_blas_test_eig","../assets/reference_lapack/TESTING/EIG","../test","stdlib_test_")


