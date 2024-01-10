# This simple script loads all Fortran77 BLAS subroutines from the reference LAPACK implementation,
# and creates a Fortran module out of them

from enum import Enum

# Create linear algebra constants module
def create_constants_module(module_name,out_folder):

    from platform import os

    INDENT          = "     "
    MAX_LINE_LENGTH = 120
    remove_headers  = True

    # Get names
    module_file = module_name + ".f90"
    module_path = os.path.join(out_folder,module_file)

    # Create file
    fid = open(module_path,"w")

    # Header
    fid.write("module {}\n".format(module_name))
    #fid.write(INDENT + "use stdlib_kinds, only: sp,dp,lk,int32,int64\n")
    fid.write(INDENT + "use iso_fortran_env, only: int32,int64\n")
    fid.write(INDENT + "use, intrinsic :: ieee_arithmetic, only: ieee_is_nan \n")
    fid.write(INDENT + "implicit none(type,external)\n")
    fid.write(INDENT + "public\n\n\n\n")

    # Temporary: to be replaced with stdlib_kinds
    fid.write(INDENT + "integer, parameter :: sp = selected_real_kind(6)\n")
    fid.write(INDENT + "integer, parameter :: dp = selected_real_kind(15)\n")
    fid.write(INDENT + "integer, parameter :: lk = kind(.true.)\n\n\n")

    # Arithmetic constants (private)
    print_lapack_constants(fid,INDENT)

    # Close module
    fid.write("\n\n\n\n\nend module {}\n".format(module_name))

    fid.close()

# Read all source files from the source folder, process them, refactor them, and put all
# subroutines/function into a module
def create_fortran_module(module_name,source_folder,out_folder,prefix,ext_functions,used_modules):

    from datetime import date
    from platform import os

    # Parameters
    today           = date.today()
    INDENT          = "     "
    MAX_LINE_LENGTH = 100
    remove_headers  = True

    # Get names
    module_file = module_name + ".f90"
    module_path = os.path.join(out_folder,module_file)

    # Get list of source files
    print("Getting list of source files...")
    source_files = []
    for file in os.listdir(source_folder):
        if (file.endswith(".f90") or file.endswith(".f") or file.endswith(".F90")) \
           and not file.startswith("la_constants") \
           and not file.startswith("la_xisnan"):
            source_files.append(file)
    source_files.sort()

    # Parse all source files
    fortran_functions = []
    for file in source_files:
        Procedures = parse_fortran_source(source_folder,file,prefix,remove_headers)
        for procedure in Procedures: fortran_functions.append(procedure)

    # Rename all procedures
    for function in fortran_functions:
        function.body,function.deps = rename_source_body(function.old_name,function.body,function.decl,\
                                                         fortran_functions,ext_functions,prefix)

    # Create file
    fid = open(module_path,"w")

    # Header
    fid.write("module {}\n".format(module_name))
    for used in used_modules:
        fid.write(INDENT + "use " + used + "\n")
    fid.write(INDENT + "implicit none(type,external)\n")
    fid.write(INDENT + "private\n\n\n\n")

    # Public interface.
    fid.write("\n\n\n" + INDENT + "public :: sp,dp,lk,int32,int64\n")
    for function in fortran_functions:
        fid.write(INDENT + "public :: " + function.new_name + "\n")

        if function.new_name=="NONAME":
            print("\n".join(function.body))
            exit(1)

    # Actual implementation
    fid.write("\n\n" + INDENT + "contains\n")

    # Write functions
    old_names,new_names = function_namelists(fortran_functions,ext_functions,prefix)
    print_function_tree(fortran_functions,old_names,fid,INDENT,MAX_LINE_LENGTH)

    # Close module
    fid.write("\n\n\nend module {}\n".format(module_name))

    fid.close()

    # Return list of all functions defined in this module, including the external ones
    return old_names

# Enum for file sections
class Section(Enum):
    HEADER = 1
    DECLARATION = 2
    EXTERNALS = 3
    BODY = 4
    END = 5

# Print LAPACK constants
def print_lapack_constants(fid,INDENT):

    real_prefix = ['s','d']
    cmpl_prefix = ['c','z']
    precision   = ['32-bit','64-bit']

    for i in range(len(real_prefix)):
       rpr = real_prefix[i]
       cpr = cmpl_prefix[i]
       rk  = rpr + "p"

       fid.write("\n" + INDENT + "! "+precision[i]+" function prefixes \n")
       fid.write(INDENT + "character,   parameter :: "+rpr+"prefix  = '"+rpr.upper()+"' \n")
       fid.write(INDENT + "character,   parameter :: "+cpr+"prefix  = '"+cpr.upper()+"' \n")

       fid.write("\n" + INDENT + "! "+precision[i]+" real constants \n")
       fid.write(INDENT + "real("+rk+"),    parameter :: "+rpr+"zero  =  0.0_"+rk+"\n")
       fid.write(INDENT + "real("+rk+"),    parameter :: "+rpr+"half  =  0.5_"+rk+"\n")
       fid.write(INDENT + "real("+rk+"),    parameter :: "+rpr+"one   =  1.0_"+rk+"\n")
       fid.write(INDENT + "real("+rk+"),    parameter :: "+rpr+"two   =  2.0_"+rk+"\n")
       fid.write(INDENT + "real("+rk+"),    parameter :: "+rpr+"three =  3.0_"+rk+"\n")
       fid.write(INDENT + "real("+rk+"),    parameter :: "+rpr+"four  =  4.0_"+rk+"\n")
       fid.write(INDENT + "real("+rk+"),    parameter :: "+rpr+"eight =  8.0_"+rk+"\n")
       fid.write(INDENT + "real("+rk+"),    parameter :: "+rpr+"ten   = 10.0_"+rk+"\n")

       fid.write("\n" + INDENT + "! "+precision[i]+" scaling constants \n")
       fid.write(INDENT + "integer,     parameter :: "     +rpr+"maxexp = maxexponent("+rpr+"zero) \n")
       fid.write(INDENT + "integer,     parameter :: "     +rpr+"minexp = minexponent("+rpr+"zero) \n")
       fid.write(INDENT + "real("+rk+"),    parameter :: "+rpr+"radix  = real(radix("+rpr+"zero),"+rk+") \n")
       fid.write(INDENT + "real("+rk+"),    parameter :: "+rpr+"ulp    = epsilon("+rpr+"zero) \n")
       fid.write(INDENT + "real("+rk+"),    parameter :: "+rpr+"eps    = "+rpr+"ulp*"+rpr+"half \n")
       fid.write(INDENT + "real("+rk+"),    parameter :: "+rpr+"safmin = "+rpr+"radix**max("+rpr+"minexp-1,1-"+rpr+"maxexp) \n")
       fid.write(INDENT + "real("+rk+"),    parameter :: "+rpr+"safmax = "+rpr+"one/"+rpr+"safmin \n")
       fid.write(INDENT + "real("+rk+"),    parameter :: "+rpr+"smlnum = "+rpr+"safmin/"+rpr+"ulp \n")
       fid.write(INDENT + "real("+rk+"),    parameter :: "+rpr+"bignum = "+rpr+"safmax*"+rpr+"ulp \n")
       fid.write(INDENT + "real("+rk+"),    parameter :: "+rpr+"rtmin  = sqrt("+rpr+"smlnum) \n")
       fid.write(INDENT + "real("+rk+"),    parameter :: "+rpr+"rtmax  = sqrt("+rpr+"bignum) \n")
       fid.write("\n" + INDENT + "! "+precision[i]+" Blue's scaling constants \n")
       fid.write(INDENT + "! ssml>=1/s and sbig==1/S with s,S as defined in https://doi.org/10.1145/355769.355771 \n")
       fid.write(INDENT + "real("+rk+"),    parameter :: "+rpr+"tsml   = "+rpr+"radix**ceiling(("+rpr+"minexp-1)*"+rpr+"half) \n")
       fid.write(INDENT + "real("+rk+"),    parameter :: "+rpr+"tbig   = "+rpr+"radix**floor(("+rpr+"maxexp-digits("+rpr+"zero)+1)*"+rpr+"half) \n")
       fid.write(INDENT + "real("+rk+"),    parameter :: "+rpr+"ssml   = "+rpr+"radix**(-floor(("+rpr+"minexp-digits("+rpr+"zero))*"+rpr+"half)) \n")
       fid.write(INDENT + "real("+rk+"),    parameter :: "+rpr+"sbig   = "+rpr+"radix**(-ceiling(("+rpr+"maxexp+digits("+rpr+"zero)-1)*"+rpr+"half)) \n")

       fid.write("\n" + INDENT + "! "+precision[i]+" complex constants \n")
       fid.write(INDENT + "complex("+rk+"), parameter :: "+cpr+"zero  = (0.0_"+rk+",0.0_"+rk+")\n")
       fid.write(INDENT + "complex("+rk+"), parameter :: "+cpr+"half  = (0.5_"+rk+",0.0_"+rk+")\n")
       fid.write(INDENT + "complex("+rk+"), parameter :: "+cpr+"one   = (1.0_"+rk+",0.0_"+rk+")\n")

    fid.write("\n\n\n" + "contains" + "\n\n\n")

# Print function tree in a dependency-suitable way
def print_function_tree(functions,fun_names,fid,INDENT,MAX_LINE_LENGTH):

    ext_funs = fun_names[len(functions):]

    # Cleanup first
    for i in range(len(functions)):
        if ext_funs is not None:
           functions[i].printed = functions[i].old_name in ext_funs
        else:
           functions[i].printed = False
        functions[i].ideps   = []

    # Get dependency indices
    for i in range(len(functions)):
        for j in range(len(functions[i].deps)):
           functions[i].ideps.append(fun_names.index(functions[i].deps[j].lower()))
        print(functions[i].ideps)

    attempt = 0
    MAXIT   = 50*len(functions)
    while attempt<MAXIT:

        attempt+=1

        for i in range(len(functions)):

            # Check deps
            if (not functions[i].printed):
                nprinted = 0
                for j in range(len(functions[i].deps)):

                    dep = functions[i].ideps[j]

                    # Dependency is an external function or the current function: do not consider
                    if dep>=len(functions) or dep==i:
                        nprinted+=1
                    elif functions[dep].printed:
                        nprinted+=1

                if nprinted==len(functions[i].deps) or attempt>=MAXIT:
                   print(str(nprinted) + "deps printed already for " + functions[i].old_name)
                   write_function_body(fid,functions[i].body,INDENT,MAX_LINE_LENGTH)
                   functions[i].printed = True


    # Final check
    not_printed = 0
    for i in range(len(functions)):
        if not functions[i].printed:
            not_printed+=1
            print(" - Function " + functions[i].old_name + "was not printed ")

    if not_printed>0:
        print("***ERROR*** there are non printed functions")
        exit(1)

# Write function body (list of lines)
def write_function_body(fid,body,INDENT,MAX_LINE_LENGTH):

    import re

    for i in range(len(body)):
       line = body[i]
       continued = False

       # Blank line
       if bool(re.match(r'^\s*!\s*$',line)):
           # If line is '!', just print a blank line
           fid.write(INDENT + "\n")
           continue

       mat = re.match(r'^\s*!', line)
       is_comment_line = bool(mat)

       # Preprocess comment line: put exclamation mark right before the first occurrence
       if is_comment_line:
          post  = line[mat.end():]
          posts = post.lstrip(' ')
          nspaces = mat.end()-mat.start()+len(post) - len(posts)
          line = (" " * nspaces) + "! " + posts

       if bool(re.match(r'^\s*!\s*$',line)):
           # If line is '!', just print a blank line
           fid.write(INDENT + "\n")
       else:

           while (len(line)>MAX_LINE_LENGTH - 2*len(INDENT)) and not is_comment_line:
              # Find last non-reserved character
              m = re.search(r'[^a-zA-Z\d\.\_\'\"\*\=\<\>\/][a-zA-Z\d\.\_\'\"\*\=\<\>\/]*$',line[:MAX_LINE_LENGTH-2])

              if m is None:
                  print(m)
                  print(line)
                  print("EEEEEEEE")
                  exit(1)

              next = line[m.start()+1:]
              end_line = "&\n" if len(next.strip())>0 else "\n"
              comment  = "continued" if continued else "non      "
              fid.write(line[:m.start()+1] + end_line)
              print(comment+" line:" + line[:m.start()+1])

              # Start with reminder (add same number of trailing spaces
              nspaces = len(line)-len(line.lstrip(' '))
              line = (" " * nspaces) + next
              print("remainder line:" + line)
              continued = True
           if len(line)>0:
               if not continued:
                   fid.write(line + "\n")
               else:
                   fid.write(INDENT + INDENT + line + "\n")

# This class represents the contents of a 1-function/1-subroutine Fortran source file parsed from BLAS/LAPACK
class Fortran_Source:
    def __init__(self):
        self.ext = ".f90"
        self.is_free_form  = True
        self.is_function   = False
        self.is_subroutine = True
        self.old_name      = "NONAME"
        self.new_name      = "NONAME"
        self.body          = []
        self.deps          = []
        self.ideps         = []
        self.decl          = []
        self.printed       = False

# Read and preprocess a Fortran line for parsing: remove comments, adjust left, and if this is a continuation
# line, read all continuation lines into it
def line_read_and_preprocess(line,is_free_form,file_name):

    import re

    processed = replace_f77_types(line,is_free_form)

    if is_free_form:
       processed = replace_la_constants(processed,file_name)

    will_continue   = bool(re.match(r".*\S+.*&\s*!*.*$", processed.rstrip()))

    if will_continue: # remove what's right of the ampersand
        print(re.sub(r'&.*\s*$','',processed).strip())
        processed = re.sub(r'&.*\s*$','',processed).strip()

    # Remove comments
    if is_free_form:
       is_comment_line = bool(re.match(r'^\s*!', processed))
       is_continuation = bool(re.match(r'^\s*&', processed))
       is_use          = bool(re.match(r'^\s*use', processed))

       # If this is a continuation line, remove all that's before the continuation character
       if is_continuation:
           processed = re.sub(r'^\s*&','',processed).strip()

    else:
       is_comment_line = bool(re.match(r'^\S\S*.*', processed))
       is_continuation = bool(re.match(r'^     [\S\&\*]', processed))

#       print("***")
#       print("processed: "+processed)
#       print(re.match(r'^     [\S\&\*]', processed))
#       print("***")

       is_use          = bool(re.match(r'^      \s*use', processed))

       # Remove continuation character
       if is_continuation:
           processed = re.sub(r'^     [\S\&\*]', '', processed).strip()

    will_continue = will_continue and not is_comment_line

    return processed,is_continuation,is_comment_line,is_use,will_continue

# Parse a line, and replace old-style Fortran datatypes and constructs with stdlib kinds
def replace_f77_types(line,is_free_form):

    import re

    INDENT = "" if is_free_form else "      "

    new_line = line.rstrip()
    new_line = new_line.replace(".LE.","<=")
    new_line = new_line.replace(".GE.",">=")
    new_line = new_line.replace(".EQ.","==")
    new_line = new_line.replace(".NE.","/=")
    new_line = new_line.replace(".LT.","<")
    new_line = new_line.replace(".GT.",">")
    new_line = new_line.replace(".AND .",".AND.")
    new_line = re.sub(r'^\s*COMPLEX\*16',INDENT+'COMPLEX(dp)',new_line)
    new_line = re.sub(r'^\s*COMPLEX ',INDENT+'COMPLEX(sp) ',new_line)
    new_line = re.sub(r'^\s*INTEGER ',INDENT+'INTEGER(int32) ',new_line)
    new_line = re.sub(r'^\s*LOGICAL ',INDENT+'LOGICAL(lk) ',new_line)
    new_line = re.sub(r'^\s*REAL ',INDENT+'REAL(sp) ',new_line)
    new_line = re.sub(r'^\s*DOUBLE PRECISION ',INDENT+'REAL(dp) ',new_line)
    new_line = re.sub(r'^\s*COMPLEX\*16,',INDENT+'COMPLEX(dp),',new_line)
    new_line = re.sub(r'^\s*COMPLEX,',INDENT+'COMPLEX(sp),',new_line)
    new_line = re.sub(r'^\s*INTEGER,',INDENT+'INTEGER(int32),',new_line)
    new_line = re.sub(r'^\s*LOGICAL,',INDENT+'LOGICAL(lk),',new_line)
    new_line = re.sub(r'^\s*REAL,',INDENT+'REAL(sp),',new_line)
    new_line = re.sub(r'^\s*DOUBLE PRECISION,',INDENT+'REAL(dp),',new_line)
    new_line = new_line.replace("E+000","_sp")
    new_line = new_line.replace("E+00","_sp")
    new_line = new_line.replace("E+01","e+01_sp")
    new_line = new_line.replace("E+0","_sp")
    new_line = new_line.replace("D+000","_dp")
    new_line = new_line.replace("D+00","_dp")
    new_line = new_line.replace("D+01","e+01_dp")
    new_line = new_line.replace("D+0","_dp")

    return new_line

def replace_la_constants(line,file_name):

    import re

    new_line = line.rstrip()

    letter = file_name[0].lower()
    if   letter=='c' or letter=='s':
        new_line = new_line.replace("_wp","_sp")
        new_line = new_line.replace("(wp)","(sp)")
        new_line = new_line.replace(" zero","szero")
        new_line = new_line.replace(" one"," sone")
        new_line = new_line.replace("(one","(sone")
        new_line = new_line.replace(" two","stwo")
        new_line = new_line.replace(" half","shalf")
    elif letter=='d' or letter=='z':
        new_line = new_line.replace("_wp","_dp")
        new_line = new_line.replace("(wp)","(dp)")
        new_line = new_line.replace(" zero","dzero")
        new_line = new_line.replace(" one"," done")
        new_line = new_line.replace("(one","(done")
        new_line = new_line.replace(" two","dtwo")
        new_line = new_line.replace(" half","dhalf")
        new_line = new_line.replace(" czero","zzero")

    return new_line

# Check if a line is a datatype line
def is_datatype_line(line):

    check_line = line.strip().lower()

    is_data_line = check_line.starts

# Check if a line is followed by external declaration
def is_externals_header(line):

    import re

    check_line = line.strip().lower()

    # Begins with a data type
    ext =    bool(re.match(r'\S\s*.. external functions ..',check_line)) \
          or bool(re.match(r'\S\s*.. external function ..',check_line)) \
          or bool(re.match(r'\S\s*external functions\s*\S*',check_line)) \
          or bool(re.match(r'\S\s*.. external subroutines ..',check_line))

    return ext

# Check if a line is a declaration line
def is_declaration_line(line):

    check_line = line.strip().lower()

    # Begins with a data type
    is_decl =    check_line.startswith("type ") \
              or check_line.startswith("type(") \
              or check_line.startswith("real ") \
              or check_line.startswith("real,") \
              or check_line.startswith("real(") \
              or check_line.startswith("real::") \
              or check_line.startswith("double precision ") \
              or check_line.startswith("double precision,") \
              or check_line.startswith("double precision(") \
              or check_line.startswith("double precision::") \
              or check_line.startswith("doubleprecision") \
              or check_line.startswith("integer ") \
              or check_line.startswith("integer,") \
              or check_line.startswith("integer(") \
              or check_line.startswith("integer::") \
              or check_line.startswith("complex ") \
              or check_line.startswith("complex,") \
              or check_line.startswith("complex(") \
              or check_line.startswith("complex::") \
              or check_line.startswith("character ") \
              or check_line.startswith("character,") \
              or check_line.startswith("character(") \
              or check_line.startswith("character*") \
              or check_line.startswith("character*(") \
              or check_line.startswith("character::") \
              or check_line.startswith("logical ") \
              or check_line.startswith("logical(") \
              or check_line.startswith("logical,") \
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
              or check_line.startswith("parameter::") \
              or check_line.startswith("implicit ")

    return is_decl

def filter_declaration_line(line):

    check_line = line.strip().lower()

    # Remove all EXTERNAL declarations
    filtered =   check_line.startswith("external ") \
              or check_line.startswith("external::") \
              or check_line.startswith("implicit ")

    return filtered

# Assemble two lists of old -> new (prefixed) function names
def function_namelists(Sources,external_funs,prefix):
    # Assemble list of old and new names
    old_names = []
    new_names = []
    for Source in Sources:
        old_names.append(Source.old_name.lower())
        new_names.append(Source.new_name.lower())

    # Add list of external functions to that of the names
    if external_funs is not None:
        for ext in external_funs:
            old_names.append(ext.lower())
            new_names.append(prefix+ext.lower())

    return old_names,new_names

# Given the list of all sources, rename all matching names in the current source body
def rename_source_body(name,lines,decl,Sources,external_funs,prefix):

    import re
    la_names = ['zero','one','two','rtmin','rtmax','safmin','safmax','sbig','ssml','tbig','tsml','half']
    la_repl = []

    initial = name[0]
    if initial=='c' or initial=='s':
        for i in range(len(la_names)):
            la_repl.append('s'+la_names[i])
    else: # z, d
        for i in range(len(la_names)):
            la_repl.append('d'+la_names[i])
        la_names.append('czero')
        la_repl .append('zzero')

    la_names.append('la_isnan')
    la_repl .append('ieee_is_nan')

    print("Renaming procedure body <"+name+">")

    body = []

    old_names,new_names = function_namelists(Sources,external_funs,prefix)

    is_found    = [False for i in range(len(new_names))]
    is_declared = [False for i in range(len(new_names))]

    la_const = False

    # First of all, map which of these names are used as declared variables. In this case
    # do not replace their names
    for i in range(len(decl)):
        old_line = decl[i].lower()
        if bool(re.search(r".+\s+[^a-zA-Z\_0-9]"+old_names[i]+r"[^a-zA-Z\_0-9].*",old_line)):
            is_declared[i] = True
        if "la_constants" in old_line:
            la_const = True

    replacement = prefix+r'\g<0>'

    whole = '\n'.join(lines)
    whole = whole.lower()

    for j in range(len(old_names)):
        if is_declared[j]: continue
        old = len(whole)
        whole = re.sub(r"\b"+old_names[j]+r"\b",replacement,whole)
        if len(whole)>old:
            print("***match***" + old_names[j])
            is_found[j] = True

    # Replace la constants
    if la_const:
        for j in range(len(la_names)):
            if is_declared[j]: continue
            old = len(whole)
            whole = re.sub(r"\b"+la_names[j]+r"\b",la_repl[j],whole)
#            if len(whole)>old:
#                print("***match***" + la_names[j])
#        print("***TEMPORARY STOP")
#        exit(1)

    body = whole.split('\n')

    # Build dependency list
    dependency_list = []
    for j in range(len(old_names)):
        if is_found[j]:
            dependency_list.append(old_names[j])

    return body,dependency_list

def parse_fortran_source(source_folder,file_name,prefix,remove_headers):

    from platform import os
    import re

    print("Parsing source file "+file_name+" ...")

    INDENT = "     "
    DEBUG  = False #file_name.lower().startswith("dlaed6")

    Procedures = []

    if file_name.endswith(".f") or file_name.endswith(".F") or file_name.endswith(".for") or file_name.endswith(".f77"):
       free_form = False
    else:
       free_form = True

    # Init empty source
    Source  = Fortran_Source()
    whereAt = Section.HEADER
    Source.is_free_form = free_form
    Source.body = []
    Source.decl = []

    # FiLoad whole file; split by lines; join concatenation lines
    with open(os.path.join(source_folder,file_name), 'r') as file:
        # Create an empty list to store the lines
        file_body = []

        # Iterate over the lines of the file
        was_continuing = False
        for line in file:
            if DEBUG: print("raw =" + line)

            # Remove the newline character at the end of the line
            line,is_continuation,is_comment,is_use,will_continue = \
               line_read_and_preprocess(line,Source.is_free_form,file_name)

            if DEBUG: print("continuation="+str(is_continuation)+" comment="+str(is_comment)+\
                            " use"+str(is_use)+" will_continue="+str(will_continue))

            # Append the line to the list, if it was not a comment line
            if was_continuing:
               if DEBUG: print("was continuing: add "+line+" to "+file_body[-1])
               file_body[-1] = file_body[-1] + line
            elif is_continuation:
               # Check if last line was a comment
               last_line,last_cont,last_comment,last_use,last_will_cont = \
               line_read_and_preprocess(file_body[-1],Source.is_free_form,file_name)

               if last_comment:
                   if DEBUG: print("last comment, add "+line+" to "+file_body[-2])
                   file_body[-2] = file_body[-2] + line
               else:
                   if DEBUG: print("last not comment, add "+line+" to "+file_body[-1])
                   file_body[-1] = file_body[-1] + line

            else:
               if DEBUG: print("new line: "+line)
               file_body.append(line)

            # Set continuation for the next line
            was_continuing = will_continue

        # Iterate over the joined lines of the file
        for line in file_body:

            if len(line.strip())<=0: continue

            # Remove the newline character at the end of the line
            line,is_continuation,is_comment,is_use,will_continue = \
            line_read_and_preprocess(line,Source.is_free_form,file_name)

            # Append the line to the list
            if is_comment:

               if DEBUG: print("Section.COMMENT " + line + " " + str(whereAt))

               # Inside an externals section: remove altogether
               if whereAt==Section.EXTERNALS:
                  if DEBUG: print("go back to declaration")
                  whereAt = Section.DECLARATION
               # Is an Externals section starting
               elif whereAt==Section.DECLARATION and is_externals_header(line):
                  if DEBUG: print("is externals header")
                  whereAt = Section.EXTERNALS
                  line = ""
               else:
                  # Just append this line, but ensure F90+ style comment
                  line = re.sub(r'^\S', '!', line)

                  if whereAt!=Section.HEADER or not remove_headers:
                     Source.body.append(INDENT + line)

            else:

               if DEBUG: print(str(whereAt) + " reads: " + line)

               # Check what section we're in
               match whereAt:
                   case Section.HEADER:
                       # Check if a declaration is starting
                       sub_found = bool('subroutine' in line.strip().lower())
                       fun_found = bool('function' in line.strip().lower())
                       name = False

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

                           if DEBUG: print("Subroutine name found: " + str(name))

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

                           if DEBUG: print("Function name found: " + str(name))

                           whereAt = Section.DECLARATION

                       # Procedure name found: add/modify header
                       if name:
                             line = line.strip()
                             if remove_headers:
                                 Source.body.append(INDENT)
                                 Source.body.append(INDENT)

                   case Section.DECLARATION:

                       # A procedure name must have been read
                       if Source.new_name=="NONAME":
                           print("INVALID PROCEDURE NAME")
                           exit(1)

                       # Check if this line still begins with a declaration
                       if is_declaration_line(line):
                           if DEBUG: print("is declaration line " + line)
                           # Filter declaration line
                           if (filter_declaration_line(line)):
                               if DEBUG: print("filter declaration line")
                               line = "";
                           else:
                               Source.decl.append(line)
                       else:
                           # Start body section
                           if DEBUG: print("start body: " + line)
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
                       if     line.strip().upper()=="END" \
                           or bool(re.match(r'^\s*END\s*SUBROUTINE.*$',line.upper())) \
                           or bool(re.match(r'^\s*END\s*FUNCTION.*$',line.upper())):
                           whereAt = Section.END
                           if Source.is_function:
                               line = "END FUNCTION " + Source.old_name.upper() + "\n"
                           elif Source.is_subroutine:
                               line = "END SUBROUTINE " + Source.old_name.upper() + "\n"

               # Append this line

               non_deleted = whereAt!=Section.HEADER or not remove_headers

               if non_deleted and not is_use:
                  Source.body.append(INDENT + line)
               else:
                  if DEBUG: print("NOT printed: " + line + " " + str(whereAt))

            # On function end
            if whereAt==Section.END:

                  # Save source
                  Procedures.append(Source)

                  # Reinitialize source
                  Source  = Fortran_Source()
                  whereAt = Section.HEADER
                  Source.is_free_form = free_form
                  Source.body = []
                  Source.decl = []

    if whereAt!=Section.END and whereAt!=Section.HEADER:
        print("WRONG SECTION REACHED!!! " + str(whereAt) + " in procedure " + Source.old_name.upper() + " file " + file_name)

    if DEBUG:
        for i in range(len(Source.body)):
           print(Source.body[i])
        exit(1)

    return Procedures



# Run script
funs = []
create_constants_module("stdlib_linalg_constants","../src")
funs = create_fortran_module("stdlib_linalg_blas",\
                             "../assets/reference_lapack/BLAS/SRC","../src",\
                             "stdlib_",\
                             funs,\
                             ["stdlib_linalg_constants"])
funs = create_fortran_module("stdlib_linalg_lapack",\
                             "../assets/reference_lapack/SRC",\
                             "../src",\
                             "stdlib_",\
                             funs,\
                             ["stdlib_linalg_constants","stdlib_linalg_blas"])
#create_fortran_module("stdlib_linalg_blas_test_eig","../assets/reference_lapack/TESTING/EIG","../test","stdlib_test_")


