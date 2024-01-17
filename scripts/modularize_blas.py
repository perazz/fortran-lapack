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
    fid.write("#if defined(_OPENMP)\n")
    fid.write(INDENT + "use omp_lib\n")
    fid.write("#endif\n")
    fid.write(INDENT + "implicit none(type,external)\n")
    fid.write(INDENT + "public\n\n\n\n")

    # Temporary: to be replaced with stdlib_kinds
    fid.write(INDENT + "integer, parameter :: sp  = selected_real_kind(6)\n")
    fid.write(INDENT + "integer, parameter :: dp  = selected_real_kind(15)\n")
    fid.write(INDENT + "integer, parameter :: lk  = kind(.true.)\n")
    fid.write(INDENT + "! Integer size support for ILP64 builds should be done here")
    fid.write(INDENT + "integer, parameter :: ilp = int32\n")
    fid.write(INDENT + "private            :: int32, int64\n\n\n")

    # Arithmetic constants (private)
    print_lapack_constants(fid,INDENT)

    # Close module
    fid.write("\n\n\n\n\nend module {}\n".format(module_name))

    fid.close()

# Patch lapack aux module interface
def patch_lapack_aux(fid,prefix,indent):

    INDENT          = "     "

    initials = ['s','d','c','z']
    datatypes = ['real(sp)','real(dp)','complex(sp)','complex(dp)']

    for i in range(len(initials)):
        fid.write(INDENT + "public :: {prf}selctg_{int}\n".format(prf=prefix,int=initials[i]))
        fid.write(INDENT + "public :: {prf}select_{int}\n".format(prf=prefix,int=initials[i]))

    fid.write("\n")


    fid.write(INDENT + "! SELCTG is a LOGICAL FUNCTION of three DOUBLE PRECISION arguments \n")
    fid.write(INDENT + "! used to select eigenvalues to sort to the top left of the Schur form. \n")
    fid.write(INDENT + "! An eigenvalue (ALPHAR(j)+ALPHAI(j))/BETA(j) is selected if SELCTG is true, i.e., \n")
    fid.write(INDENT + "abstract interface \n")
    for i in range(len(initials)):
        if (i<=1):
            fid.write(INDENT + "   logical(lk) function {prf}selctg_{int}(alphar,alphai,beta) \n".format(prf=prefix,int=initials[i]))
            fid.write(INDENT + "       import sp,dp,lk \n")
            fid.write(INDENT + "       implicit none \n")
            fid.write(INDENT + "       {}, intent(in) :: alphar,alphai,beta \n".format(datatypes[i]))
            fid.write(INDENT + "   end function {prf}selctg_{int} \n".format(prf=prefix,int=initials[i]))
            fid.write(INDENT + "   logical(lk) function {prf}select_{int}(alphar,alphai) \n".format(prf=prefix,int=initials[i]))
            fid.write(INDENT + "       import sp,dp,lk \n")
            fid.write(INDENT + "       implicit none \n")
            fid.write(INDENT + "       {}, intent(in) :: alphar,alphai \n".format(datatypes[i]))
            fid.write(INDENT + "   end function {prf}select_{int} \n".format(prf=prefix,int=initials[i]))
        else:
            fid.write(INDENT + "   logical(lk) function {prf}selctg_{int}(alpha,beta) \n".format(prf=prefix,int=initials[i]))
            fid.write(INDENT + "       import sp,dp,lk \n")
            fid.write(INDENT + "       implicit none \n")
            fid.write(INDENT + "       {}, intent(in) :: alpha,beta \n".format(datatypes[i]))
            fid.write(INDENT + "   end function {prf}selctg_{int} \n".format(prf=prefix,int=initials[i]))
            fid.write(INDENT + "   logical(lk) function {prf}select_{int}(alpha) \n".format(prf=prefix,int=initials[i]))
            fid.write(INDENT + "       import sp,dp,lk \n")
            fid.write(INDENT + "       implicit none \n")
            fid.write(INDENT + "       {}, intent(in) :: alpha \n".format(datatypes[i]))
            fid.write(INDENT + "   end function {prf}select_{int} \n".format(prf=prefix,int=initials[i]))


    fid.write(INDENT + "end interface \n\n")


# Read all source files from the source folder, process them, refactor them, and put all
# subroutines/function into a module
def create_fortran_module(module_name,source_folder,out_folder,prefix,ext_functions,used_modules, \
                          split_by_initial):

    from datetime import date
    from platform import os


    # Splitting by initials
    if split_by_initial:
        initials = ['aux','s','d','c','z']
        modules = []
        for i in range(len(initials)):
            modules.append(module_name+"_"+initials[i])
    else:
        modules = [module_name]
        initials = ['']

    # Parameters
    today           = date.today()
    INDENT          = "     "
    MAX_LINE_LENGTH = 100
    remove_headers  = True

    # Get list of source files
    print("Getting list of source files...")
    source_files = []
    for file in os.listdir(source_folder):
        if (file.endswith(".f90") or file.endswith(".f") or file.endswith(".F90") or file.endswith(".F")) \
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

    # Create modules
    for m in range(len(modules)):
        this_module = module_name
        if len(initials[m])>0: this_module = this_module + "_" + initials[m]
        module_file = this_module + ".f90"
        module_path = os.path.join(out_folder,module_file)
        fid = open(module_path,"w")

        # Header
        fid.write("module {}\n".format(this_module))
        for used in used_modules:
            fid.write(INDENT + "use " + used + "\n")

        # Add top modules in the hierarchy
        for n in range(m):
            fid.write(INDENT + "use " + module_name + "_" + initials[n] + "\n")

        fid.write(INDENT + "implicit none(type,external)\n")
        fid.write(INDENT + "private\n\n\n\n")

        # Public interface.
        fid.write("\n\n\n" + INDENT + "public :: sp,dp,lk,ilp\n")
        for function in fortran_functions:
            if function_in_module(initials[m],function.old_name):
                fid.write(INDENT + "public :: " + function.new_name + "\n")

            if function.new_name=="NONAME":
                print("\n".join(function.body))
                exit(1)

        # Patch AUX
        if module_name+"_"+initials[m]=='stdlib_linalg_lapack_aux':
            patch_lapack_aux(fid,prefix,INDENT)

        # Actual implementation
        fid.write("\n\n" + INDENT + "contains\n")

        # Write functions
        old_names,new_names = function_namelists(fortran_functions,ext_functions,prefix)
        print_function_tree(fortran_functions,old_names,fid,INDENT,MAX_LINE_LENGTH,initials[m])

        # Close module
        fid.write("\n\n\nend module {}\n".format(this_module))
        fid.close()

    # Write wrapper module
    if split_by_initial:
        module_file = module_name + ".f90"
        module_path = os.path.join(out_folder,module_file)

        fid = open(module_path,"w")

        # Header
        fid.write("module {}\n".format(module_name))
        for used in used_modules:
            fid.write(INDENT + "use " + used + "\n")

        for i in initials:
            fid.write(INDENT + "use {mname}_{minit}\n".format(mname=module_name,minit=i))
        fid.write(INDENT + "implicit none(type,external)\n")
        fid.write(INDENT + "public\n")
        # Close module
        fid.write("\n\n\nend module {}\n".format(module_name))
        fid.close()


    # Return list of all functions defined in this module, including the external ones
    return old_names

def function_module_initial(function_name):
   initials = ['aux','c','s','d','z']

   oname = function_name.lower().strip()

   for i in range(len(initials)):
       initial = initials[i]
       if len(initial)<1:
           return 'a'
       elif    oname.endswith("amax") \
            or oname.endswith("abs") \
            or oname.endswith("abs1") :
           return 'a'
       elif initial in ['c','s','d','z'] and oname[0]==initial[0].lower():
               return oname[0]

   # No matches
   return 'a'

def function_in_module(initial,function_name):

   oname = function_name.lower().strip()

   if len(initial)<1:
       in_module = True
   elif    oname.endswith("amax") \
        or oname.endswith("abs") \
        or oname.endswith("roundup_lwork") \
        or oname.endswith("abs1") :
       in_module = initial[0].lower() == 'a'
   # PATCH: exclude functions
   # - with names ending in *x or *_extended, as they require external subroutines
   # which are not provided by the Fortran implementation
   elif ((len(oname)>6 and oname.endswith('x')) or \
         (len(oname)>6 and (oname.endswith('2') \
                            and not oname.endswith('ladiv2')   \
                            and not oname.endswith('geqrt2')   \
                            and not oname.endswith('getrf2')   \
                            and not oname.endswith('getrfnp2') \
                            and not oname.endswith('tplqt2')   \
                            and not oname.endswith('sytrs2')   \
                            and not oname.endswith('hetrs2')   \
                            and not oname.endswith('orbdb2')   \
                            and not oname.endswith('unbdb2')   \
                            and not oname.endswith('tpqrt2')   \
                            and not oname.endswith('potrf2'))) or \
         (len(oname)>6 and (oname.endswith('3') \
                            and not oname.endswith('ladiv3')   \
                            and not oname.endswith('geqrt3')   \
                            and not oname.endswith('tpqrt3')   \
                            and not oname.endswith('gelqt3')   \
                            and not oname.endswith('sytrs_3')  \
                            and not oname.endswith('hetrs_3')  \
                            and not oname.endswith('orbdb3')   \
                            and not oname.endswith('unbdb3')   \
                            and not oname.endswith('trevc3'))) \
          or oname.endswith('extended') \
          or oname.endswith('ssytri2') \
          or oname.endswith('_2stage') \
          or oname.endswith('ssysv_rk')):
       in_module = False
   elif initial[0].lower() in ['c','s','d','z'] :
       in_module = oname[0]==initial[0].lower()
   else:
       in_module = not (oname[0] in ['c','s','d','z'])
   return in_module

# Enum for file sections
class Section(Enum):
    HEADER = 1
    HEADER_DESCR = 2
    DECLARATION = 3
    EXTERNALS = 4
    BODY = 5
    END = 6


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
def print_function_tree(functions,fun_names,fid,INDENT,MAX_LINE_LENGTH,initial):

    ext_funs = fun_names[len(functions):]

    # Cleanup first. Mark all functions that will go into another module as "printed",
    # so they won't be checked as dependencies
    for i in range(len(functions)):
        functions[i].ideps   = []
        if ext_funs is not None:
           functions[i].printed = functions[i].old_name in ext_funs
        else:
           functions[i].printed = False

        functions[i].printed = functions[i].printed or \
                               not function_in_module(initial,functions[i].old_name)

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
                   write_function_body(fid,functions[i].header," " * header_indentation(functions[i].body),MAX_LINE_LENGTH,False)
                   write_function_body(fid,functions[i].body,INDENT,MAX_LINE_LENGTH,True)
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

# Check if line is directive
def is_directive_line(line):

   ls = line.lstrip()

   is_dir = line.startswith("#") or \
            line.startswith("!$OMP") or \
            line.startswith("!$omp")

   return is_dir

# Return indentation for the header
def header_indentation(body):
    nspaces = 1
    # Seek SUBROUTINE or FUNCTION line
    for i in range(len(body)):
        l = body[i].lower()
        if 'subroutine' in l or 'function' in l:
            nspaces = heading_spaces(body[i])
            print("INDENT FOUND: "+str(nspaces)+" FROM <"+l+">")
            break

    return nspaces

# Return number of spaces at the beginning of a string
def heading_spaces(line):
    posts = line.lstrip(' ')
    nspaces = len(line)-len(posts)
    return nspaces

# Adjust variable declaration
def adjust_variable_declaration(line,datatype):

    import re

    declarations = [r'integer\(\S+\)', \
                    r'character\*\(\S+\)', \
                    r'real\(\S+\)', \
                    r'complex\(\S+\)', \
                    r'logical\(\S+\)', \
                    r'character', \
                    r'intrinsic', \
                    r'external']

    ll = line.lower()

    for i in range(len(declarations)):
        m  = re.match(r'^\s*' + declarations[i] + r'\s+\w',ll)
        if m:
            variable = line[m.end()-1:].lstrip()
            line = line[:m.end()-2].rstrip() + " :: " + variable

            # Patch function argument
            if i==4 and variable.lower().strip()=='selctg':
                print(line)
                nspaces = len(line)-len(line.lstrip(' '))
                line = nspaces*" " + "procedure(stdlib_selctg_"+datatype[0]+") :: selctg"
                print(line)

            if i==4 and variable.lower().strip()=='select':
                print(line)
                nspaces = len(line)-len(line.lstrip(' '))
                line = nspaces*" " + "procedure(stdlib_select_"+datatype[0]+") :: select"
                print(line)

            return line

    return line

# Write function body (list of lines)
def write_function_body(fid,body,INDENT,MAX_LINE_LENGTH,adjust_comments):

    import re

    fid.write("\n")

    header = True

    for i in range(len(body)):
       line = body[i]
       continued = False

       # Blank line
       if line.strip()=="":
           # Do not print
           if header: continue
       else:
           header = False

       # Blank comment line
       if bool(re.match(r'^\s*!\s*\.{0,2}\s*$',line)):
           # If line is '!', just skip it
           continue

       # Patches
       find = [r'larfg\( n, a, a\( 1, min\( 2, n \) \), lda, t \)$']
       repl = [r'larfg( n, a(1,1), a( 1, min( 2, n ) ), lda, t(1,1) )']

       for j in range(len(find)):
           line = re.sub(find[j],repl[j],line)


       mat = re.match(r'^\s*!', line)
       is_comment_line = bool(mat)
       is_directive = is_directive_line(line)

       # Preprocess comment line: put exclamation mark right before the first occurrence
       if is_comment_line and not is_directive:
          if adjust_comments:
              # Put exclamation mark at the location where the first nonspace character was
              post  = line[mat.end():]
              posts = post.lstrip(' ')
              nspaces = mat.end()-mat.start()+len(post) - len(posts)
              line = (" " * nspaces) + "! " + posts
          else:
              # Just add indent
              line = INDENT + line.lstrip(' ')

       if is_directive:
           fid.write(line+"\n")
       elif bool(re.match(r'^\s*!\s*$',line)):
           # If line is '!', just print a blank line
           fid.write(INDENT + "\n")
       else:

           while (len(line)>MAX_LINE_LENGTH - 2*len(INDENT)) and not is_comment_line:

              shift = 0

              # Find last non-reserved character
              m = re.search(r'[^a-zA-Z\d\.\_\'\"\*\=\<\>\/][a-zA-Z\d\.\_\'\"\*\=\<\>\/]*$',line[:MAX_LINE_LENGTH-2])

              if m is None:
                  print(m)
                  print(line)
                  print("EEEEEEEE")
                  exit(1)

              # PATCH :: Check that we're not splitting a string between quotes, aka ' &\n'
              if re.search(r'\'\s+$',line[:m.start()+1]): shift = -2

              next = line[m.start()+1+shift:]

              end_line = "&\n" if len(next.strip())>0 else "\n"
              comment  = "continued" if continued else "non      "
              fid.write(line[:m.start()+1+shift] + end_line)
              print(comment+" line:" + line[:m.start()+1+shift])

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

class Fortran_Line:
    def __init(self):
        self.string        = ""
        self.continuation  = False
        self.comment       = False
        self.use           = False
        self.will_continue = False
        self.directive     = False

# Read and preprocess a Fortran line for parsing: remove comments, adjust left, and if this is a continuation
# line, read all continuation lines into it
def line_read_and_preprocess(line,is_free_form,file_name):

    import re

    processed = replace_f77_types(line,is_free_form)

    if is_free_form:
       processed = replace_la_constants(processed,file_name)

    # Check if this is a directive
    is_dir = is_directive_line(processed)

    will_continue   = bool(re.match(r".*\S+.*&\s*!*.*$", processed.rstrip()))

    if will_continue and not is_dir: # remove what's right of the ampersand
        processed = re.sub(r'&.*\s*$','',processed).strip()

    # Remove comments
    if is_free_form:
       is_comment_line = bool(re.match(r'^\s*!', processed))
       is_continuation = bool(re.match(r'^\s*&', processed))
       is_use          = bool(re.match(r'^\s*use', processed))

       # If this is a continuation line, remove all that's before the continuation character
       if is_continuation and not is_dir:
           processed = re.sub(r'^\s*&','',processed).strip()

    else:
       is_comment_line = bool(re.match(r'^\S\S*.*', processed))
       is_continuation = bool(re.match(r'^     [\S\&\*]', processed))


       is_use          = bool(re.match(r'^      \s*use', processed))

       # Remove continuation character
       if is_continuation and not is_dir:
           processed = re.sub(r'^     [\S\&\*]', '', processed).strip()

    will_continue = will_continue and not is_comment_line

    Line = Fortran_Line()
    Line.string = processed
    Line.continuation = is_continuation
    Line.comment = is_comment_line
    Line.use = is_use
    Line.will_continue = will_continue
    Line.directive = is_dir

    return Line

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
    new_line = re.sub(r'^\s*INTEGER ',INDENT+'INTEGER(ilp) ',new_line)
    new_line = re.sub(r'^\s*LOGICAL ',INDENT+'LOGICAL(lk) ',new_line)
    new_line = re.sub(r'^\s*REAL ',INDENT+'REAL(sp) ',new_line)
    new_line = re.sub(r'^\s*DOUBLE PRECISION ',INDENT+'REAL(dp) ',new_line)
    new_line = re.sub(r'^\s*COMPLEX\*16,',INDENT+'COMPLEX(dp),',new_line)
    new_line = re.sub(r'^\s*COMPLEX,',INDENT+'COMPLEX(sp),',new_line)
    new_line = re.sub(r'^\s*INTEGER,',INDENT+'INTEGER(ilp),',new_line)
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
    new_line = new_line.replace("0.0d0","0.0_dp")

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
          or bool(re.match(r'\S\s*from blas',check_line)) \
          or bool(re.match(r'\S\s*from lapack',check_line)) \
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
    whole_decl = '\n'.join(decl).lower()
    for j in range(len(old_names)):
        if bool(re.search(r"\b"+old_names[j]+r"\b",whole_decl)): is_declared[i] = True
        if "la_constants" in whole_decl: la_const = True

    replacement = prefix+r'\g<0>'


    whole = '\n'.join(lines).lower()
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

    body = whole.split('\n')

    # Restore directive lines cases
    for j in range(len(body)):
       if is_directive_line(body[j]): body[j] = lines[j]

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

    initial = 'a'

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
    Source.header = []
    open_loops = []
    loop_starts = []

    # FiLoad whole file; split by lines; join concatenation lines
    with open(os.path.join(source_folder,file_name), 'r') as file:
        # Create an empty list to store the lines
        file_body = []

        # Iterate over the lines of the file
        was_continuing = False
        for line in file:
            if DEBUG: print("raw =" + line)

            # Remove the newline character at the end of the line
            Line = line_read_and_preprocess(line,Source.is_free_form,file_name)

            if DEBUG: print("continuation="+str(Line.continuation)+" comment="+str(Line.comment)+\
                            " use"+str(Line.use)+" will_continue="+str(Line.will_continue))

            # This is a directive: append as-is
            if Line.directive:
                file_body.append(Line.string)
            # Append the line to the list, if it was not a comment line
            elif was_continuing:
               if DEBUG: print("was continuing: add "+Line.string+" to "+file_body[-1])
               file_body[-1] = file_body[-1] + Line.string
            elif Line.continuation:
               # Check if last line was a comment
               Last_Line = line_read_and_preprocess(file_body[-1],Source.is_free_form,file_name)

               if Last_Line.comment:
                   if DEBUG: print("last comment, add "+Line.string+" to "+file_body[-2])
                   file_body[-2] = file_body[-2] + Line.string
               else:
                   if DEBUG: print("last not comment, add "+Line.string+" to "+file_body[-1])
                   file_body[-1] = file_body[-1] + Line.string

            else:
               if DEBUG: print("new line: "+Line.string)
               file_body.append(Line.string)

            # Set continuation for the next line
            was_continuing = Line.will_continue

        # Iterate over the joined lines of the file
        for line in file_body:

            if len(line.strip())<=0: continue

            # Remove the newline character at the end of the line
            Line = line_read_and_preprocess(line,Source.is_free_form,file_name)

            # Append the line to the list
            if Line.directive:

               # Directives: apend as-is
               Source.body.append(line)

            elif Line.comment:

               if DEBUG: print("Section.COMMENT " + line + " " + str(whereAt))

               # Inside an externals section: remove altogether
               if whereAt==Section.EXTERNALS:
                  if DEBUG: print("go back to declaration")
                  whereAt = Section.DECLARATION
               # Is an Externals section starting
               elif whereAt==Section.DECLARATION and is_externals_header(line):
                  if DEBUG: print("is externals header")
                  whereAt = Section.EXTERNALS
                  continue
               # We are parsing the header description line
               elif whereAt==Section.HEADER_DESCR:
                  if DEBUG: print("is description header")

                  # Go back to standard header
                  if '\endverbatim' in line:
                      whereAt=Section.HEADER
                  else:
                      # Strip comment sign and add to header
                      line = re.sub(r'\*\>',"",line)
                      line = re.sub(r'^\*',"",line)
                      line = re.sub(r'\!\>',"",line)
                      line = re.sub(r'\\verbatim',"",line)
                      line = re.sub(r'=============',"",line)
                      if len(line)>0: Source.header.append("! " + line.strip())
                      if DEBUG: print(Source.header[-1])


               else:
                  # Just append this line, but ensure F90+ style comment
                  line = re.sub(r'^\S', '!', line)

                  if whereAt!=Section.HEADER or not remove_headers:
                     Source.body.append(INDENT + line)
                  elif remove_headers:
                     # Remove header, only keep description
                     if '\par Purpose:' in line:
                        whereAt = Section.HEADER_DESCR

               # Empty comment line? skip
               ls = line.strip()
               if ls=='!' or ls=='! ..': continue

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
                               initial = function_module_initial(Source.old_name)

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
                               initial = function_module_initial(Source.old_name)

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

                               line = adjust_variable_declaration(line,initial)

                               Source.decl.append(line)
                       else:
                           # Start body section
                           if DEBUG: print("start body: " + line)
                           whereAt = Section.BODY

                   case Section.EXTERNALS:

                       # Empty line after a comment? skip it

                       # Check if this line still begins with a declaration
                       if is_declaration_line(line):
                           # Delete altoghether
                           line = "";
                       else:
                           # Go back to declaration section
                           whereAt = Section.DECLARATION

                   case Section.BODY:

                       # Labelled do loop? save
                       m_loop = re.match(r'^\s*do\s+[0-9]+\s+',line.lower())
                       if m_loop:
                           # Extract label
                           numbers = re.findall(r'\d+',line)
                           nspaces = len(line) - len(line.lstrip(' '))
                           line    = (" "*nspaces) + "loop_" + str(numbers[0]) + ": do " + line[m_loop.end():]

                           open_loops.append(numbers[0])
                           loop_starts.append(nspaces)

                       # Go to inside labelled loop
                       m_goto = re.search(r'go\s*to\s+\d+$',line.lower())
                       if bool(m_goto):
                           # Extract label
                           numbers = re.findall(r'\d+$',line.lower())

                           # This "go to" matches one of the current loops
                           if len(open_loops)>0:
                               if numbers[0] in open_loops:
                                   nspaces = len(line) - len(line.lstrip(' '))
                                   left = line[:m_goto.start()]
                                   if len(left.strip())<=0:
                                       line = (nspaces*" ") + "cycle loop_" + str(numbers[0])
                                   else:
                                       line = line[:m_goto.start()] + "cycle loop_" + str(numbers[0])

                                   if DEBUG: print(line)

                       # End of labelled loop
                       if re.match(r'^\s+\d+\s+continue',line.lower()):
                           # Extract label
                           numbers = re.findall(r'\d+',line.lower())

                           # This "continue" matches loop
                           if len(open_loops)>0:
                               if open_loops[-1]==numbers[0]:
                                   loop_ID = open_loops.pop()
                                   nspaces = loop_starts.pop()
                                   line    = (" "*nspaces) + "end do loop_" + str(loop_ID)

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
               non_use     = whereAt!=Section.DECLARATION or not Line.use

               if non_deleted and non_use:
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
                  Source.body   = []
                  Source.decl   = []
                  Source.header = []
                  open_loops    = []
                  loop_starts   = []

    if whereAt!=Section.END and whereAt!=Section.HEADER:
        print("WRONG SECTION REACHED!!! " + str(whereAt) + " in procedure " + Source.old_name.upper() + " file " + file_name)

    if DEBUG:
        for i in range(len(Source.body)):
           print(Source.body[i])
        exit(1)

    return Procedures


# Copy files into a temporary folder
import shutil
import os
import glob

os.makedirs('../assets/lapack_sources',exist_ok=True)


for file in glob.glob(r'../assets/reference_lapack/SRC/*.f*') + glob.glob(r'../assets/reference_lapack/SRC/*.F*'):
    print(file)
    shutil.copy(file, '../assets/lapack_sources/')

shutil.copyfile('../assets/reference_lapack/INSTALL/slamch.f', '../assets/lapack_sources/slamch.f')
shutil.copyfile('../assets/reference_lapack/INSTALL/dlamch.f', '../assets/lapack_sources/dlamch.f')
shutil.copyfile('../assets/reference_lapack/INSTALL/sroundup_lwork.f', '../assets/lapack_sources/sroundup_lwork.f')
shutil.copyfile('../assets/reference_lapack/INSTALL/droundup_lwork.f', '../assets/lapack_sources/droundup_lwork.f')

# Run script
funs = []
create_constants_module("stdlib_linalg_constants","../src")
funs = create_fortran_module("stdlib_linalg_blas",\
                             "../assets/reference_lapack/BLAS/SRC","../src",\
                             "stdlib_",\
                             funs,\
                             ["stdlib_linalg_constants"],True)
funs = create_fortran_module("stdlib_linalg_lapack",\
                             "../assets/lapack_sources",\
                             "../src",\
                             "stdlib_",\
                             funs,\
                             ["stdlib_linalg_constants","stdlib_linalg_blas"],True)
#create_fortran_module("stdlib_linalg_blas_test_eig","../assets/reference_lapack/TESTING/EIG","../test","stdlib_test_")




