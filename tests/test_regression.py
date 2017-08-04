#!/usr/bin/env python
# MARTINIZE
# A simple, versatile tool for coarse-graining molecular systems
# Copyright (C) 2017  Tsjerk A. Wassenaar and contributors
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA  02110-1301, USA.

"""
Regression tests for martinize.

This test suite runs the insane command line with various set of arguments, and
assess that the results correspond to the result obtained with previous
versions.

Notice that these tests do not assess that the results are correct. Instead,
they assess that changes do not affect the behaviour of the program.

If ran as a script, this generate the reference files expected by the tests. If
ran usinf pytest or nosetest, this executes insane with a series of arguments
and compares the output to the reference.
"""

from __future__ import print_function

import contextlib
import functools
import glob
import mock
import os
import random
import shutil
import shlex
import subprocess
import sys
import tempfile
import textwrap
import importlib

import testfixtures
from nose.tools import assert_equal, assert_raises

import utils

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

try:
    from itertools import izip_longest as zip_longest
except ImportError:
    from itertools import zip_longest


PROGRAM = 'martinize'

CLI = importlib.import_module(".".join([PROGRAM, "cli"]))

HERE = os.path.abspath(os.path.dirname(__file__))
EXECUTABLE = utils.which(PROGRAM)
DATA_DIR = os.path.join(HERE, 'data')
INPUT_DIR = os.path.join(HERE, 'data', 'inputs')
RANDSEED = '42'
SEED_ENV = 'INSANE_SEED'

PDB_LIST = ('1ubq', '3csy', '2qwo', '1a8g', '2oar')  #, '1cag')
FF_LIST = ('martini21', 'martini21p',
           'martini22', 'martini22p',
           'elnedyn', 'elnedyn22', 'elnedyn22p')

# The arguments to test insane with are listed here. The tuple is used both to
# generate the references, and to run the tests.
# To add a test case, add the arguments to test in the tuple.
SIMPLE_TEST_CASES = [
    ('-f {}.pdb'.format(pdb), pdb) for pdb in PDB_LIST
]
SIMPLE_TEST_CASES.extend([
    # Examples from the martini tutorial
    # <http://cgmartini.nl/index.php/tutorials-general-introduction-gmx5/proteins-gmx5>
    ('-f 1ubq.pdb -o system-vaccum.top -x 1UBQ-CG.pdb '
     '-dssp dssp -p backbone -ff martini22', '1ubq'),
    ('-f 1ubq.pdb -o system-vaccum.top -x 1UBQ-CG.pdb '
     '-ss chainA.ss -p backbone -ff martini22', '1ubq-ss'),
    ('-f 1a8g.pdb -o system-vaccum.top -x 1A8G-CG.pdb '
     '-dssp dssp -p backbone -ff martini22', '1a8g'),
    ('-f 1a8g.pdb -o system-vaccum.top -x 1A8G-CG.pdb '
     '-dssp dssp -p backbone -ff martini22 '
     '-elastic -ef 500 -el 0.5 -eu 0.9 -ea 0 -ep 0', '1a8g', '1a8g-elastic'),
    ('-f 1a8g.pdb -o system-vaccum.top -x 1UBQ-CG.pdb '
     '-dssp dssp -p backbone -ff elnedyn22', '1a8g'),
])
SIMPLE_TEST_CASES.extend([
    # Examples taken from Djurre's tests
    # <https://github.com/cgmartini/martinize.py/blob/master/test/test.sh>
    ('-f 1ubq.pdb -o 1UBQ_cg.top -x 1UBQ_cg.pdb '
     '-ss ~EEEEEETTS~EEEEE~~TTSBHHHHHHHHHHHH~~~GGGEEEEETTEE~~TTSBTGGGT~~TT~EEEEEE~~S~~',
     '1ubq', '1ubq-inline-ss'),
    ('-f 2oar.pdb -o 2OAR_cg.top -x 2OAR_cg.pdb '
     '-sep -nt -p All -pf 500 -dssp dssp -ff martini22', '2oar'),
    ('-f 1cag.pdb -o 1CAG_cg.top -x 1CAG_cg.pdb -collagen -ff martini22', '1cag'),
    ('-f 3sjm.pdb -o 3SJM_cg.top -x 3SJM_cg.pdb -collagen -ff martini22dna', '3sjm'),
])
SIMPLE_TEST_CASES.extend([
    ('-f 1l35.pdb -o 1L35_cg.top -x 1l35_cg.pdb '
     '-cys auto -name lysozyme -dssp dssp -ed -ff {}'.format(ff), '1l35')
    for ff in FF_LIST
])
SIMPLE_TEST_CASES.extend([
    ('-f 1a8g.pdb -merge A,B', '1a8g'),
    ('-f 2oar.pdb -merge A,B,C -merge D,E', '2oar'),
])
SIMPLE_TEST_CASES.extend([
    ('-f {}.pdb -ff {}'.format(pdb, ff), pdb)
    for pdf in PDB_LIST
    for ff in FF_LIST
])
SIMPLE_TEST_CASES.extend([
    ('-f {}.pdb -nmap nmap.idx'.format(pdb), pdb) for pdb in PDB_LIST
])
SIMPLE_TEST_CASES.extend([
    ('-f {}.pdb -n index.idx'.format(pdb), pdb) for pdb in PDB_LIST
])


def _arguments_as_list(arguments):
    """
    Return the arguments as a list as expected by subprocess.Popen.

    The arguments can be provided as a string that will be spitted to a list.
    They can also be provided as a list, then the list will be returned
    untouched.
    """
    try:
        arguments_list = shlex.split(arguments)
    except ValueError:
        arguments_list = arguments
    return arguments_list


def _split_case(case):
    """
    Get the arguments and the input directory from a test case.
    """
    if len(case) == 3:
        case_args, input_dir, alias = case
        if input_dir is not None:
            input_dir = os.path.join(INPUT_DIR, input_dir)
    elif len(case) == 2:
        case_args, input_dir = case
        input_dir = os.path.join(INPUT_DIR, input_dir)
        alias = case_args
    else:
        case_args = case
        input_dir = None
        alias = case_args
    return case_args, input_dir, alias


def _reference_path(arguments, alias=None):
    """
    Get the path to the reference files for the simple test cases.
    """
    arg_list = _arguments_as_list(arguments)
    simple_case_ref_data = os.path.join(DATA_DIR, 'simple_case')
    base_name = arguments if alias is None else alias
    return os.path.join(simple_case_ref_data, base_name)


def _run_external(arguments):
    command = [EXECUTABLE] + arguments
    env = {SEED_ENV: INSANE_SEED} if SEED_ENV else {}
    process = subprocess.Popen(command,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               env=env)
    out, err = process.communicate()
    log = ''
    return out, err, log, process.returncode


def _run_internal(arguments):
    if SEED_ENV:
        os.environ[SEED_ENV] = RANDSEED
    random.seed(RANDSEED)
    command = [EXECUTABLE] + arguments
    out = StringIO()
    err = StringIO()
    with mock.patch('random.random', return_value=0.1):
        with utils._redirect_out_and_err(out, err):
            with testfixtures.LogCapture() as log:
                returncode = CLI.main(command)
    out = out.getvalue()
    err = err.getvalue()
    log = str(log)
    return out, err, log, returncode


def run_program(arguments, input_directory=None, runner=_run_internal):
    """
    Run program with the given arguments

    Insane is run in a copy of `input_directory`.
    """
    # Copy the content of the input directory in the current directory if an
    # input directory is provided.
    if input_directory is not None:
        for path in glob.glob(os.path.join(input_directory, '*')):
            if os.path.isdir(path):
                shutil.copytree(path, '.')
            else:
                shutil.copy2(path, '.')
    out, err, log, returncode = runner(arguments)
    print("** {} exited with return code {}.".format(PROGRAM.capitalize(), returncode))
    if returncode:
        print(err)
    return out, err, log, returncode


def run_and_compare(arguments, input_dir, ref_dir, runner):
    """
    Run program and compare its output against a reference
    """
    # Create the command as a list for subprocess.Popen.
    # The arguments can be pass to the current function as a string or as a
    # list of arguments. If they are passed as a string, they need to be
    # converted to a list.
    arguments = _arguments_as_list(arguments)

    ref_stdout = os.path.join(ref_dir, 'stdout')
    ref_stderr = os.path.join(ref_dir, 'stderr')
    ref_log = os.path.join(ref_dir, 'testlog')

    # We want program to run in a temporary directory. This allows to keep the
    # file system clean, and it avoids mixing output of different tests.
    with utils.tempdir():
        out, err, log, returncode = run_program(arguments, input_dir, runner=runner)
        assert not returncode
        utils.compare(utils.ContextStringIO(out), ref_stdout)
        utils.compare(utils.ContextStringIO(err), ref_stderr)
        utils.compare(utils.ContextStringIO(log), ref_log)
        utils.compare_directories('./', ref_dir,
                                  ignore=('stderr', 'stdout', 'testlog'))


def _test_simple_cases():
    """
    This function generates test functions for nosetests. These test functions
    execute insane with the argument listed in SIMPLE_TEST_CASES.
    """
    for case in SIMPLE_TEST_CASES:
        case_args, input_dir, alias = _split_case(case)
        ref_gro, ref_top, ref_stdout, ref_stderr = _reference_path(case_args, alias)
        # The test generator could yield run and compare directly. Bt, then,
        # the verbose display of nosetests gets crowded with the very long
        # names of the reference file, that are very redundant. Using a partial
        # function allows to have only the arguments for insane displayed.
        _test_case = functools.partial(
            run_and_compare,
            ref_gro=ref_gro,
            ref_top=ref_top,
            ref_stdout=ref_stdout,
            ref_stderr=ref_stderr,
            runner=_run_internal)
        _test_case.__doc__ = ' '.join([EXECUTABLE, case_args])
        yield (_test_case, case_args, input_dir)


def test_simple_cases_internal():
    """
    This function generates test functions for nosetests. These test functions
    calls insane's main function with the argument listed in SIMPLE_TEST_CASES.
    """
    for case in SIMPLE_TEST_CASES:
        case_args, input_dir, alias = _split_case(case)
        ref_dir = _reference_path(case_args, alias)
        # The test generator could yield run and compare directly. Bt, then,
        # the verbose display of nosetests gets crowded with the very long
        # names of the reference file, that are very redundant. Using a partial
        # function allows to have only the arguments for insane displayed.
        _test_case = functools.partial(
            run_and_compare,
            ref_dir=ref_dir,
            runner=_run_internal)
        yield (_test_case, case_args, input_dir)


class TestGroTester(object):
    """
    Test if the comparison of GRO file catches the differences.
    """
    ref_gro_content = """\
    INSANE! Membrane UpperLeaflet>POPC=1 LowerLeaflet>POPC=1
    4
        1POPC   NC3    1   2.111  14.647  11.951
        1POPC   PO4    2   2.177  14.644  11.651
        1POPC   GL1    3   2.128  14.642  11.351
        1POPC   GL2    4   1.961  14.651  11.351
    10 10 10"""

    def test_equal(self):
        """
        Make sure that identical files do not fail.
        """
        with utils.tempdir():
            with open('ref.gro', 'w') as outfile:
                print(textwrap.dedent(self.ref_gro_content),
                      file=outfile, end='')
            utils.assert_gro_equal('ref.gro', 'ref.gro')

    def test_diff_x(self):
        """
        Make sure that error in coordinates is caught.
        """
        gro_content = """\
        INSANE! Membrane UpperLeaflet>POPC=1 LowerLeaflet>POPC=1
        4
            1POPC   NC3    1   2.111  14.647  11.951
            1POPC   PO4    2   2.177  14.644  11.651
            1POPC   GL1    3   2.128  14.642  11.353  # Is not within tolerance
            1POPC   GL2    4   1.961  14.651  11.351
        10 10 10"""

        with utils.tempdir():
            with open('ref.gro', 'w') as outfile:
                print(textwrap.dedent(self.ref_gro_content),
                      file=outfile, end='')
            with open('content.gro', 'w') as outfile:
                print(textwrap.dedent(gro_content), file=outfile, end='')
            assert_raises(AssertionError, utils.assert_gro_equal,
                          'content.gro', 'ref.gro')

    def test_diff_in_tolerance(self):
        """
        Make sure that small errors in coordinates are not caught.
        """
        gro_content = """\
        INSANE! Membrane UpperLeaflet>POPC=1 LowerLeaflet>POPC=1
        4
            1POPC   NC3    1   2.111  14.647  11.951
            1POPC   PO4    2   2.177  14.644  11.651
            1POPC   GL1    3   2.128  14.642  11.352  # Is within tolerance
            1POPC   GL2    4   1.961  14.651  11.351
        10 10 10"""

        with utils.tempdir():
            with open('ref.gro', 'w') as outfile:
                print(textwrap.dedent(self.ref_gro_content),
                      file=outfile, end='')
            with open('content.gro', 'w') as outfile:
                print(textwrap.dedent(gro_content), file=outfile, end='')
            utils.assert_gro_equal('content.gro', 'ref.gro')

    def test_diff_natoms(self):
        """
        Make sure that differences in number of atom is caught.
        """
        gro_content = """\
        INSANE! Membrane UpperLeaflet>POPC=1 LowerLeaflet>POPC=1
        6
            1POPC   NC3    1   2.111  14.647  11.951
            1POPC   PO4    2   2.177  14.644  11.651
            1POPC   GL1    3   2.128  14.642  11.351
            1POPC   GL2    4   1.961  14.651  11.351
            1POPC   C1A    5   2.125  14.651  11.051
            1POPC   D2A    6   2.134  14.602  10.751
        10 10 10"""

        with utils.tempdir():
            with open('ref.gro', 'w') as outfile:
                print(textwrap.dedent(self.ref_gro_content),
                      file=outfile, end='')
            with open('content.gro', 'w') as outfile:
                print(textwrap.dedent(gro_content), file=outfile, end='')
            assert_raises(AssertionError, utils.assert_gro_equal,
                          'content.gro', 'ref.gro')

    def test_diff_title(self):
        """
        Make sure that a different title is caught.
        """
        gro_content = """\
        A different title
        4
            1POPC   NC3    1   2.111  14.647  11.951
            1POPC   PO4    2   2.177  14.644  11.651
            1POPC   GL1    3   2.128  14.642  11.351
            1POPC   GL2    4   1.961  14.651  11.351
        10 10 10"""

        with utils.tempdir():
            with open('ref.gro', 'w') as outfile:
                print(textwrap.dedent(self.ref_gro_content),
                      file=outfile, end='')
            with open('content.gro', 'w') as outfile:
                print(textwrap.dedent(gro_content), file=outfile, end='')
            assert_raises(AssertionError, utils.assert_gro_equal,
                          'content.gro', 'ref.gro')

    def test_diff_box(self):
        """
        Make sure that a different box is caught.
        """
        gro_content = """\
        INSANE! Membrane UpperLeaflet>POPC=1 LowerLeaflet>POPC=1
        4
            1POPC   NC3    1   2.111  14.647  11.951
            1POPC   PO4    2   2.177  14.644  11.651
            1POPC   GL1    3   2.128  14.642  11.351
            1POPC   GL2    4   1.961  14.651  11.351
        10 9.9 10 9.08 4 54"""

        with utils.tempdir():
            with open('ref.gro', 'w') as outfile:
                print(textwrap.dedent(self.ref_gro_content),
                      file=outfile, end='')
            with open('content.gro', 'w') as outfile:
                print(textwrap.dedent(gro_content), file=outfile, end='')
            assert_raises(AssertionError, utils.assert_gro_equal,
                          'content.gro', 'ref.gro')

    def test_diff_field(self):
        """
        Make sure that a difference in a field is caught.
        """
        gro_content = """\
        INSANE! Membrane UpperLeaflet>POPC=1 LowerLeaflet>POPC=1
        4
            1POPC   NC3    1   2.111  14.647  11.951
            1DIFF   PO4    2   2.177  14.644  11.651
            1POPC   GL1    3   2.128  14.642  11.351
            1POPC   GL2    4   1.961  14.651  11.351
        10 10 10"""

        with utils.tempdir():
            with open('ref.gro', 'w') as outfile:
                print(textwrap.dedent(self.ref_gro_content),
                      file=outfile, end='')
            with open('content.gro', 'w') as outfile:
                print(textwrap.dedent(gro_content), file=outfile, end='')
            assert_raises(AssertionError, utils.assert_gro_equal,
                          'content.gro', 'ref.gro')


def generate_simple_case_references():
    """
    Run program to generate reference files for the simple regression tests.
    """
    for case in SIMPLE_TEST_CASES:
        case_args, input_dir, alias = _split_case(case)
        arguments = _arguments_as_list(case_args)
        ref_dir = _reference_path(case_args, alias)
        ref_stdout = os.path.join(ref_dir, 'stdout')
        ref_stderr = os.path.join(ref_dir, 'stderr')
        ref_log = os.path.join(ref_dir, 'testlog')

        if os.path.exists(ref_dir):
            shutil.rmtree(ref_dir)
        os.mkdir(ref_dir)

        with utils.in_directory(ref_dir):
            print(PROGRAM + ' ' + ' '.join(arguments))
            out, err, log, _ = run_program(arguments, input_dir)
            with open(ref_stdout, 'w') as outfile:
                for line in out:
                    print(line, file=outfile, end='')
            with open(ref_stderr, 'w') as outfile:
                for line in err:
                    print(line, file=outfile, end='')
            with open(ref_log, 'w') as outfile:
                print(log, end='', file=outfile)



def clean_simple_case_references():
    """
    Delete reference files for the simple tests if they are not in use anymore.
    """
    simple_test_cases = [_split_case(case)[2] for case in SIMPLE_TEST_CASES]
    simple_case_ref_data = os.path.join(DATA_DIR, 'simple_case')
    for path in glob.glob(os.path.join(simple_case_ref_data, '*')):
        base_name = os.path.basename(os.path.splitext(path)[0])
        if base_name not in simple_test_cases:
            print(path)
            os.remove(path)


def main():
    """
    Command line entry point.
    """
    help_ = """
Generate or clean the reference files for program's regression tests.

{0} gen: generate the files
{0} clean: clean the unused files

nosetests -v: run the tests
""".format(sys.argv[0])
    commands = {'gen': generate_simple_case_references,
                'clean': clean_simple_case_references}
    if len(sys.argv) != 2:
        print(help_, file=sys.stderr)
        sys.exit(1)
    try:
        commands[sys.argv[1]]()
    #except KeyError:
    #    print("Unrecognized keyword '{}'.".format(sys.argv[1]))
    #    print(help_, file=sys.stderr)
    #    sys.exit(1)
    finally:
        pass


if __name__ == '__main__':
    main()
