import invoke
import os
import glob

# invoke build-pybind11 cppfile ---> build cppfile for binding
# python pybind11_example.py ---> run example program


@invoke.task
def clean(c):
    """Remove any built objects"""
    file_patterns = [
        "*.o",
        "*.so",
        "*.obj",
        "*.dll",
        "*.exp",
        "*.lib",
        "*.pyd",
        "cffi_example*",
        "cython_wrapper.cpp",
    ]
    for file_pattern in file_patterns:
        for file in glob.glob(file_pattern):
            os.remove(file)


@invoke.task()
def build_cppcode(cpp_input, cpp_output):
    """Build a shared library for C++ code"""
    invoke.run(
        "g++ -O3 -Wall -Werror -shared -std=c++11 -fPIC {0} "
        "-o {1} ".format(cpp_input, cpp_output)
    )
    print("* Complete")


def compile_python_module(cpp_name, extension_name):
    invoke.run(
        "g++ -O3 -shared -std=c++11 -fPIC "
        "`python3 -m pybind11 --includes` "
        "-I /usr/include/python3.7 -I .  "
        "{0} "
        "-o {1}`python3.8-config --extension-suffix` "
        " -Wl,-rpath,.".format(cpp_name, extension_name)
        # "-L. -lcppmult -Wl,-rpath,.".format(cpp_name, extension_name)
    )


@invoke.task()
def build_pybind11(c, cppfile):
    """Build the pybind11 wrapper library"""
    # provide cpp file name without extension
    compile_python_module(f"{cppfile}.cpp", cppfile)
    print("* Complete")


@invoke.task()
def test_pybind11(c, test_file):
    """Run the script to test PyBind11"""
    invoke.run("python3 {0}.py".format(test_file), pty=True)


@invoke.task(
    clean,
    build_pybind11,
    test_pybind11,
)
def all(c):
    """Build and run all tests"""
    pass
