#!/bin/sh

echo "#include <string.h>\n\n"         > $2
echo "int api_str_to_int( char * string, int * value) {" >> $2
grep "[ID]PARM_[A-Z_]*[ ]*=[ ]*[0-9]+" $1 | sed -e 's/\([ID]PARM_[A-Z_]*\)[ ]*=[ ]*\([0-9]*\).*/  if(!strcmp("\1", string)) { *value = \2\; return 0\;}/' >> $2
grep "API_[A-Z_]*[ ]*=[ ]*[0-9]*"      $1 | sed -e 's/\(API_[A-Z_]*\)[ ]*=[ ]*\([0-9]*\).*/  if(!strcmp("\1", string)) { *value = \2\; return 0\;}/'      >> $2
echo "  return 1;" >> $2
echo "}" >> $2

