"""
Wrappers
========

 @file wrappers/__init__.py

 PaStiX wrapper generators module intialization

 @copyright 2017-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.4.0
 @author Mathieu Faverge
 @author Selmane Lebdaoui
 @author Tony Delarue
 @date 2024-07-05

"""

# exclude inline functions from the interface
exclude_list = [ "spmIntSort1Asc1", "spmIntSort2Asc1",
                 "spmIntSort2Asc2", "spmIntMSortIntAsc",
                 "spmScalMatrix", "spmScalVector",
                 "orderDraw", "orderSupernodes",
                 "pastixOrderCompute", "pastixOrderApplyLevelOrder",
                 "pastixOrderExpand", "pastixOrderAmalgamate",
                 "pastixOrderAddIsolate", "pastixOrderFindSupernodes" ]

# translation_table with names of auxiliary variables
return_variables_dict = {
    "int":            ("info"),
    "double":         ("dval"),
    "float":          ("fval"),
    "pastix_int_t":   ("ival"),
    "spm_int_t":      ("ival"),
    "char":           ("retval"),
    "void":           ("retval"),
    "pastix_order_t": ("order"),
    "pastix_rhs_t":   ("rhs"),
    "spmatrix_t":     ("spmo"),
}

# global list used to determine derived types
derived_types = [ 'spmatrix_t', 'spm_int_t', 'pastix_int_t', 'pastix_data_t', 'pastix_order_t', 'pastix_rhs_t', 'MPI_Comm' ]

# name arrays which will be translated to assumed-size arrays, e.g. pA(*)
arrays_names_2D = ["pA", "pB", "pC", "pAB", "pQ", "pX", "pAs", "A", "B", "C", "Bl", "Bg", "X", "X0", "S" ]
arrays_names_1D = ["colptr", "rowptr", "loc2glob", "dofs", "row", "values",
                   "iparm", "dparm", "bindtab", "perm", "invp", "schur_list",
                   "rang", "tree", "x0", "x", "y", "b", "list" ]

def polish_file(whole_file):
    """Preprocessing and cleaning of the header file.
       Do not change the order of the regular expressions !
       Works with a long string."""

    clean_file = whole_file

    # borrowed from cfwrapper.py
    # Remove C comments:
    clean_file = re.sub(r"(?s)/\*.*?\*/", "", clean_file)
    clean_file = re.sub(r"//.*", "", clean_file)
    # Remove BEGIN/END_C_DECLS statement:
    clean_file = re.sub("BEGIN_C_DECLS", "", clean_file)
    clean_file = re.sub("END_C_DECLS", "", clean_file)
    # Remove C directives (multilines then monoline):
    clean_file = re.sub(r"(?m)^#(.*[\\][\n])+.*?$", "", clean_file)
    clean_file = re.sub("(?m)^#.*$", "", clean_file)
    clean_file = re.sub("(?m)#.*", "", clean_file)
    # Remove TABs and overnumerous spaces:
    clean_file = clean_file.replace("\t", " ")
    clean_file = re.sub("[ ]{2,}", " ", clean_file)
    # Remove extern C statement:
    clean_file = re.sub("(?m)^(extern).*$", "", clean_file)
    # Remove empty lines:
    clean_file = re.sub(r"(?m)^\n$", "", clean_file)
    # Merge structs
    clean_file = re.sub(r"(?m)$", "", clean_file)

    # Merge string into single line
    clean_file = re.sub(r"\n", "", clean_file)

    # Split the line based on ";" and "}"
    clean_file = re.sub(r";", "\n", clean_file)
    clean_file = re.sub(r"}", "}\n", clean_file)

    clean_file = re.sub(r"\bSPM_Comm", "MPI_Comm", clean_file)
    clean_file = re.sub(r"\bPASTIX_Comm", "MPI_Comm", clean_file)

    return clean_file

def preprocess_list(initial_list):
    """Preprocessing and cleaning of the header file.
       Works with a list of strings.
       Produces a new list in which each function, enum or struct
       corresponds to a single item."""

    # merge braces
    list1 = []
    merged_line = ""
    nopen = 0
    inStruct = False
    for line in initial_list:

        if (line.find("struct") > -1):
            inStruct = True

        if (inStruct):
            split_character = ","
        else:
            split_character = ""

        nopen += line.count("{") - line.count("}")
        merged_line += line + split_character

        if (nopen <= 0):
            list1.append(merged_line)
            merged_line = ""
            isOpen   = False
            inStruct = False
            nopen = 0

    # merge structs
    list2 = []
    merged_line = ""
    for line in list1:

        merged_line += line

        if (line.find("struct") == -1):
            list2.append(merged_line)
            merged_line = ""

    # clean orphan braces
    list3 = []
    for line in list2:

        if (line.strip() != "}"):
            list3.append(line)

    #print '\n\n'.join(list3)

    return list3

def parse_triple( string ):
    """Parse string of
       type (*)name
       into triple of [type, pointer, name, const]"""

    string = string.strip()

    # Split name from type
    m = re.search(r"^(.*[\s\*])([^\*\s]+)$", string )
    if m == None:
        print("Error: Cannot detect type for ", string)

    namestr = m.group(2).strip()
    pointer = ""
    const   = 0
    typestr = m.group(1).strip()

    # Remove useless const on parameters
    typestr = re.sub( r"const$", "", typestr ).strip()

    if "const" in string:
        const = 1
        typestr = re.sub( r"const", "", typestr )

    # Remove spaces between *
    typestr = re.sub( r"\*\s*\*", "**", typestr )

    # Double pointer case
    if typestr[-2:] == "**":
        typestr = typestr[:-2].strip()
        pointer = "**"

    if typestr[-1:] == "*":
        typestr = typestr[:-1].strip()
        pointer = "*"

    return [typestr, pointer, namestr, const]

def parse_arg( string ):
    """Parse string of
       type (*)name
       into triple of [type, pointer, name, const]"""

    string = string.strip()

    # Split name from type
    m = re.search(r"^(.*[\s\*])([^\*\s]+)$", string )
    if m == None:
        if re.search(r"^\s*void\s*$", string ):
            return None
            namestr = ""
            typestr = "void"
        else:
            print("Error: Cannot detect type for ", string)
            return None
    else:
        namestr = m.group(2).strip()
        typestr = m.group(1).strip()

    arg = { 'name'    : namestr,
            'const'   : False,
            'pointer' : 0,
            'type'    : None }

    # Remove useless const on parameters
    typestr = re.sub( r"const$", "", typestr ).strip()

    if "const" in string:
        arg['const'] = True
        typestr = re.sub( r"const\s*", "", typestr )

    # Remove spaces between *
    typestr = re.sub( r"\*\s*\*", "**", typestr )

    # Double pointer case
    if typestr[-2:] == "**":
        typestr = typestr[:-2].strip()
        arg['pointer'] = 2

    if typestr[-1:] == "*":
        typestr = typestr[:-1].strip()
        arg['pointer'] = 1

    arg['type'] = typestr
    return arg

def parse_enums(preprocessed_list):
    """Each enum will be parsed into a list of its arguments."""

    enum_list = []
    for proto in preprocessed_list:

        # extract the part of the function from the prototype
        fun_parts = proto.split("{")

        split_fun = fun_parts[0].strip().split()
        if len(split_fun) == 0:
            continue

        if ((split_fun[0] == "enum") or
            (split_fun[0] == "typedef" and split_fun[1] == "enum")):


            if split_fun[0] == "enum":
                enumname = split_fun[1]
            else:
                enumname = split_fun[2]
            enumname = re.sub(r"_e$", "", enumname)
            enumname = re.sub(r"^pastix_", "", enumname)
            enumname = re.sub(r"^spm_", "", enumname)

            args_string = fun_parts[1];
            args_string = re.sub(r"}", "", args_string)
            args_list = args_string.split(",")
            params = [];
            for args in args_list:
                args = args.strip();
                if (args != ""):
                    values = args.split("=")

                    name = values[0].strip()
                    if (len(values) > 1):
                        value = values[1].strip()
                    else:
                        if (len(params) > 0):
                            value = params[len(params)-1][1] + 1
                        else:
                            value = 0

                    params.append([name, value])

            enum_list.append([enumname, params])

    return enum_list


def parse_structs(preprocessed_list):
    """Each struct will be parsed into a list of its arguments."""

    struct_list = []
    for proto in preprocessed_list:

        # extract the part of the function from the prototype
        fun_parts = proto.split("{")

        if (fun_parts[0].find("struct") > -1) and (len(fun_parts) > 1):
            args_string = fun_parts[1]
            parts = args_string.split("}")
            args_string = parts[0].strip()
            args_string = re.sub(r"volatile", "", args_string)
            if (len(parts) > 1):
                name_string = parts[1]
                name_string = re.sub(r"(?m),", "", name_string)
                name_string = name_string.strip()
            else:
                print("Error: Cannot detect name for ", proto)
                name_string = "name_not_recognized"

            args_list = args_string.split(",")
            params = [];
            params.append(["struct","",name_string])
            for arg in args_list:
                if (not (arg == "" or arg == " ")):
                    params.append(parse_triple(arg))

            struct_list.append(params)
            derived_types.append(name_string)

    # reorder the list so that only defined types are exported
    goAgain = True
    while (goAgain):
        goAgain = False
        for istruct in range(0,len(struct_list)-1):
            struct = struct_list[istruct]
            for j in range(1,len(struct)-1):
                type_name = struct_list[istruct][j][0]

                if (type_name in wrappers.derived_types):

                    # try to find the name in the registered types
                    definedEarlier = False
                    for jstruct in range(0,istruct):
                        struct2 = struct_list[jstruct]
                        that_name = struct2[0][2]
                        if (that_name == type_name):
                            definedEarlier = True

                    # if not found, try to find it behind
                    if (not definedEarlier):
                        definedLater = False
                        for jstruct in range(istruct+1,len(struct_list)-1):
                            struct2 = struct_list[jstruct]
                            that_name = struct2[0][2]
                            if (that_name == type_name):
                                index = jstruct
                                definedLater = True

                        # swap the entries
                        if (definedLater):
                            print("Swapping " + struct_list[istruct][0][2] + " and " + struct_list[index][0][2])
                            tmp = struct_list[index]
                            struct_list[index] = struct_list[istruct]
                            struct_list[istruct] = tmp
                            goAgain = True
                        else:
                            print("Error: Cannot find a derived type " + type_name + " in imported structs.")

    return struct_list


def parse_prototypes( preprocessed_list ):
    """Each prototype will be parsed into a list of its arguments."""

    function_list = []
    for proto in preprocessed_list:

        # Not a function decalration
        if (proto.find("(") == -1):
            continue

        fun_parts = proto.split("(")
        if len(fun_parts) < 2:
            print( "Error: The function does not have aguments", fun_def )
            continue

        fun_def  = str.strip( fun_parts[0] )
        fun_args = str.strip( fun_parts[1] )

        #
        # Extract the return type and the function name
        #
        m = re.search( r"^(.*[\s\*])([^\*\s]+)$", fun_def )
        if m == None:
            print( "Error: Cannot detect function name and return type for:  ", fun_def )

        function = {}
        function['name']    = m.group(2).strip()
        function['rettype'] = m.group(1).strip()

        # Filter and clean the return type
        if function['rettype'].find("inline") != -1:
            continue

        if function['rettype'].find("static") != -1:
            print( "Found a static non inlined function. Can't generate a wrapper for: " + function['name'] )
            continue

        # Skip excluded functions
        if function['name'] in exclude_list:
            continue

        function['rettype'] = parse_arg( function['rettype'] + " " + function['name'] )
        # Is it a function or a subroutine ?
        function['is_function'] = False
        if (function['rettype']['pointer'] != 0) or (function['rettype']['type'] != "void"):
            function['is_function'] = True

        #
        # Process the arguments
        #
        fun_args = fun_args.split(")")[0]
        fun_args_list = fun_args.split(",")

        # generate argument list
        args_list = []
        # append arguments
        for arg in fun_args_list:
            arg = arg.strip();
            if arg == "":
                continue
            cleanarg = parse_arg( arg )
            if cleanarg != None:
                args_list.append( cleanarg )

        function['args'] = args_list

        # add it only if there is no duplicity with previous one
        is_function_already_present = False
        for fun in function_list:
            if function['name'] == fun['name']:
                is_function_already_present = True
                break

        if not is_function_already_present:
            function_list.append( function )

    return function_list

def gen_enum_copy( input_enum ):
    output_enum = list()
    output_enum.append( input_enum[0] )
    output_enum.append( list() )

    for e in input_enum[1]:
        output_enum[1].append( e[:] )

    return output_enum

__all__ = [ 'exclude_list', 'return_variables_dict', 'derived_types', 'arrays_names_1D', 'arrays_names_2D',
            'polish_file', 'preprocess_list', 'parse_triple', 'parse_enums', 'parse_structs',
            'parse_prototypes', 'gen_enum_copy' ]

from .wrap_python  import *
from .wrap_fortran import *
from .wrap_julia   import *

