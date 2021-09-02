"""
 @file gen_parm_2_csv.py

 Declaration of the pastix_param2csv routine in parse_option.c.

 @copyright 2021-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.2.0
 @author Tony Delarue
 @date 2021-04-07

"""

documentation = '''/**
 *******************************************************************************
 *
 * @brief Dump the iparm an dparm parameters in the CSV file.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The main data structure.
 *
 * @param[inout] csv
 *          The csv file that will contain the dumped datas.
 *
 *******************************************************************************/
'''

def genParm2csv( iparms, dparms ) :
    result  = documentation
    result += '''void
pastix_param2csv( const pastix_data_t *pastix_data,
                        FILE          *csv )
{
    pastix_int_t *iparm = pastix_data->iparm;
    double       *dparm = pastix_data->dparm;

'''
    iline = r'    fprintf( csv, "%s,%ld\n", '
    dline = r'    fprintf( csv, "%s,%e\n",  '
    eline = r'    fprintf( csv, "%s,%s\n",  '
    for groups in iparms :
        currMaxSize = max( list( map(lambda x : len(x['name']), groups["subgroup"]) ) )
        for iparm in groups["subgroup"] :
            name   = iparm["name"]
            spaces = " " * ( currMaxSize - len(name) )

            if "enum" in iparm :
                result += eline + "\"" + name + "\", " + spaces + " pastix_"+ iparm["enum"] +"_getstr(iparm[" + name.upper() + "]) );\n"
            else :
                result += iline + "\"" + name + "\","  + spaces + " (long)iparm[" + name.upper() + "] );\n"
        result += "\n"

    currMaxSize = max( list( map(lambda x : len(x['name']), dparms) ) )
    for dparm in dparms :
        name   = dparm["name"]
        spaces = " " * ( currMaxSize - len(name) )
        result += dline + "\"" + name + "\"," + spaces + " dparm[" + name.upper() + "] );\n"

    result += "}\n"

    return result
