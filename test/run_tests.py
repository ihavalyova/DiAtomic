import sys, os
import optparse
import unittest

def getoptionparser():

    parser = optparse.OptionParser()

    parser.add_option("-q", "--quiet",
                      action="store_const", const=0, dest="verbose", default=1,
                      help="minimal output")
    parser.add_option("-v", "--verbose",
                      action="store_const", const=2, dest="verbose", default=1,
                      help="verbose output")
    parser.add_option("-i", "--include", type="string",
                      action="append",  dest="include", default=[],
                      help="include tests matching PATTERN", metavar="PATTERN")
    parser.add_option("-e", "--exclude", type="string",
                      action="append", dest="exclude", default=[],
                      help="exclude tests matching PATTERN", metavar="PATTERN")
    parser.add_option("-f", "--failfast",
                      action="store_true", dest="failfast", default=False,
                      help="stop on first failure")
    parser.add_option("-c", "--catch",
                      action="store_true", dest="catchbreak", default=False,
                      help="catch Control-C and display results")

    return parser

def import_package(options, pkgname):

    package = __import__(pkgname)

    return package

def setup_unittest(options):

    from unittest import TestSuite
    try:
        from unittest.runner import _WritelnDecorator
    except ImportError:
        from unittest import _WritelnDecorator

    writeln_orig = _WritelnDecorator.writeln

    def writeln(self, message=''):
        try: 
            self.stream.flush()
        except: 
            pass
        writeln_orig(self, message)

        try: 
            self.stream.flush()
        except: 
            pass

    _WritelnDecorator.writeln = writeln

def getbuilddir():

    from distutils.util import get_platform
    s = os.path.join("build", "lib.%s-%.3s" % (get_platform(), sys.version))
    if hasattr(sys, 'gettotalrefcount'): s += '-pydebug'
    return s

def setup_python(options):

    rootdir = os.path.dirname(os.path.dirname(__file__))
    builddir = os.path.join(rootdir, getbuilddir())

def getpythoninfo():

    x, y = sys.version_info[:2]
    return ("Python %d.%d (%s)" % (x, y, sys.executable))

def writeln(message='', endl='\n'):

    sys.stderr.flush()
    sys.stderr.write(message+endl)
    sys.stderr.flush()

def print_banner(options, package):

    fmt = "[%d@%s] %s"
    if options.verbose:
        writeln(getpythoninfo())
        #writeln(fmt % (getlibraryinfo()))
        #writeln(fmt % (getpackageinfo(package)))

def load_tests(options, args):

    # Find tests
    import re, glob
    testsuitedir = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, testsuitedir)
    pattern = 'test_*.py'
    wildcard = os.path.join(testsuitedir, pattern)
    testfiles = glob.glob(wildcard)
    include = exclude = None

    if options.include:
        include = re.compile('|'.join(options.include)).search
    if options.exclude:
        exclude = re.compile('|'.join(options.exclude)).search
    testnames = []
    for testfile in testfiles:
        filename = os.path.basename(testfile)
        testname = os.path.splitext(filename)[0]
        if ((exclude and exclude(testname)) or
            (include and not include(testname))):
            continue
        testnames.append(testname)
    testnames.sort()

    # Load tests and populate suite
    testloader = unittest.TestLoader()
    testsuite = unittest.TestSuite()
    for testname in testnames:
        module = __import__(testname)
        for arg in args:
            try:
                cases = testloader.loadTestsFromNames((arg,), module)
            except AttributeError:
                continue
            testsuite.addTests(cases)
        if not args:
            cases = testloader.loadTestsFromModule(module)
            testsuite.addTests(cases)
    return testsuite

def run_tests(options, testsuite, runner=None):

    if runner is None:
        runner = unittest.TextTestRunner()
        runner.verbosity = options.verbose
        runner.failfast = options.failfast

    if options.catchbreak:
        unittest.installHandler()

    result = runner.run(testsuite)
    return result.wasSuccessful()

def main(args=None):

    pkgname = 'DiAtomic'

    parser = getoptionparser()
    options, args = parser.parse_args(args)
    setup_python(options)
    setup_unittest(options)
    package = import_package(options, pkgname)
    print_banner(options, package)
    testsuite = load_tests(options, args)
    success = run_tests(options, testsuite)
    
    return not success

if __name__ == '__main__':
    import sys
    sys.dont_write_bytecode = True
    sys.exit(main())
