# Copyright (C) 2011 Pierre de Buyl

# This file is part of pyMPCD

# pyMPCD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# pyMPCD is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with pyMPCD.  If not, see <http://www.gnu.org/licenses/>.

from numpy.distutils.core import setup, Extension
import os
import sys
import subprocess

ISRELEASED = False
MAJOR = 0
MINOR = 1
MICRO = 0
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

# Return the git revision as a string
def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout = subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip()
    except OSError:
        GIT_REVISION = "Unknwn"

    return GIT_REVISION

def write_version_py(filename='pyMPCD/version.py'):
    cnt = """# Copyright (C) 2011 Pierre de Buyl

# This file is part of pyMPCD

# pyMPCD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# pyMPCD is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with pyMPCD.  If not, see <http://www.gnu.org/licenses/>.

# THIS FILE IS GENERATED FROM pyMPCD SETUP.PY
short_version='%(version)s'
version='%(version)s'
full_version='%(full_version)s'
git_revision='%(git_revision)s'
release=%(isrelease)s

if not release:
    version = full_version
"""
    FULL_VERSION = VERSION
    if not ISRELEASED:
        FULL_VERSION += '.dev'
        if os.path.exists('.git'):
            GIT_REVISION = git_version()
        elif os.path.exists(filename):
            # must be a source distribution, use existing version file
            from pyMPCD.version import git_revision as GIT_REVISION
        else:
            GIT_REVISION = "Unknwn"

        FULL_VERSION += GIT_REVISION[:6]

    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION,
                       'full_version' : FULL_VERSION,
                       'git_revision' : GIT_REVISION,
                       'isrelease': str(ISRELEASED)})
    finally:
        a.close()
    return FULL_VERSION

if sys.version_info[0] < 3:
    import __builtin__ as builtins
else:
    import builtins

# From numpy setup
builtins.__PYMPCD_SETUP__ = True

full_version = write_version_py()

m1 = Extension('pyMPCD.MPCD_f', sources = [ 'pyMPCD/MPCD_f.f90' ] ,
               extra_compile_args=['-DF2PY_REPORT_ON_ARRAY_COPY=1'])

setup(name="pyMPCD", version=full_version, 
      description="pyMPCD - Multiparticle Collision Dynamics", 
      author="Pierre de Buyl", author_email="john@doe.com", 
      maintainer="Pierre de Buyl", 
      maintainer_email="john@doe.com", 
      license="GPLv3", url="http://homepages.ulb.ac.be/~pdebuyl/pympcd/", 
      packages=["pyMPCD"],
      ext_modules = [m1])

