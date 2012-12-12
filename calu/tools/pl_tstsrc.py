#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys

G_F = "Fortran source"
G_C = "C source"
G_H = "C header"

G_F_CONT_CHAR = "$" # continuation character in Fortran
G_F_COMM_CHAR = "*" # comment character in Fortran
G_F_MAX_LLEN  = 80  # maximum number of characters per line

class SourceFile:
  def __init__(self, fname):
    self.fname = fname
    f = open(fname)
    self.data = f.read()
    f.close()
    if fname[-2:] == ".f":
      self.ftpe = G_F
    elif fname[-2:] == ".c":
      self.ftpe = G_C
    elif fname[-2:] == ".h":
      self.ftpe = G_H

    self.wcount = 0

  def is_c_file(self):
    if self.ftpe == G_C or self.ftpe == G_H:
      return True
    return False

  def is_fortran_file(self):
    if self.ftpe == G_F:
      return True
    return False

  def check_last_line(self):
    if self.data[-1] != "\n":
      self.warn("No EOL at EOF")

  def check_CR(self):
    if self.data.find("\r") >= 0:
      self.warn("CR in file")

  def check_TAB(self):
    if self.data.find("\t") >= 0:
      self.warn("TAB in file")

  def check_wspace(self):
    for lno0, line in enumerate(self.data.split("\n")):
      if len(line) > 0 and line[-1].isspace():
        self.warn("White space at the end of line %d" % (lno0+1))

  def check_f_comments(self):
    for lno0, line in enumerate(self.data.split("\n")):
      if len(line) > 0 and line[0] != " " and line[0] != G_F_COMM_CHAR:
        self.warn("Wrong comment character '%c' in line %d" % (line[0], lno0+1))

  def check_f_continuation(self):
    for lno0, line in enumerate(self.data.split("\n")):
      if len(line) > 5 and line[0] == " " and line[5] != " " and line[5] != G_F_CONT_CHAR:
        self.warn("Wrong continuation character '%c' in line %d" % (line[5], lno0+1))

  def check_f_llength(self):
    for lno0, line in enumerate(self.data.split("\n")):
      if len(line) > G_F_MAX_LLEN:
        self.warn("Line %d too long: %d." % (lno0+1, len(line)))

  def check_f_do(self):
    for lno0, line in enumerate(self.data.split("\n")):
      if len(line) > 0 and line[0] == " ":
        tkns = line.lower().split()
        if "do" == tkns[0]:
          try:
            label = int(tkns[1])
          except:
            self.warn("Wrong DO syntax in line %d: %s" % (lno0+1, line))

  def check_f_badkeys(self):
    for k in ("ENDDO",):
      if k in self.data:
        self.warn("Wrong keyword: %s" % k)

  def check_c_comment(self):
    if self.data.find("//"):
      self.warn("Found // in the source code.")

  def warn(self, s):
    self.wcount += 1
    print s

def check_file(fname):
  print "Checking", fname
  sf = SourceFile(fname)
  sf.check_last_line()
  sf.check_CR()
  sf.check_TAB()
  sf.check_wspace()
  if sf.is_fortran_file():
    sf.check_f_llength()
    sf.check_f_comments()
    sf.check_f_continuation()
    sf.check_f_do()
    sf.check_f_badkeys()
  if sf.is_c_file():
    sf.check_c_comment()

def main(argv):
  for fname in argv[1:]:
    check_file(fname)

  return 0

if "__main__" == __name__:
  sys.exit(main(sys.argv))