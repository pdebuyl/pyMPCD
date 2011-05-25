# -*- coding: utf-8 -*-

# Copyright (C) 2011 Pierre de Buyl

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 only.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
suppress_docstrings.py

This scripts removes all docstrings from a Python file. It is intended to be 
used as a filter for doxygen documentation.

Usage: python suppress_docstrings.py file.py
The result is output on stdout.

Any file whose name does not end by .py is given back as is.
"""

from sys import argv

f = file(argv[1], 'r')

if (argv[1][-3:]!='.py'):
    for l in f.readlines():
        print l,
    exit()

import re

# find a class or function declaration
my_def = re.compile(r'(?:\s*)(def|class)(?:\s+)(\w+)\(')
# find a one line docstring
d_oneline = re.compile(r'^(\s*)(\"\"\"|\'\'\').*(\"\"\"|\'\'\')(\s*)$')
# find the beginning of a multiline docstring
db = re.compile(r'^(\s*)(\"\"\"|\'\'\')')
# find the end of a multiline docstring
de = re.compile(r'.*(\"\"\"|\'\'\')(\s*)$')

# beginning of file
BOF = -1
# in file
INF = 0
# in definition
IND = 1
# in multiline docstring
INM = 2

state = -1
for l in f.readlines():
    if (state==BOF):
        no = d_oneline.match(l)
        if no:
            state=INF
        else:
            n = db.match(l)
            if n:
                state=INM
            else:
                print l,
            continue
    m = my_def.match(l)
    if m:
        state = IND
        print l, 
    else:
        if (state==IND):
            no = d_oneline.match(l)
            if no:
                state=INF
            else:
                n = db.match(l)
                if n:
                    state=INM
                else:
                    print l,
        elif (state==INM):
            n = de.match(l)
            if n: 
                state=INF
        else:
            print l,

f.close()
