import os
import numpy
import tempfile
import subprocess
from setuptools import setup, Extension, find_packages
from distutils import log
from distutils.command.install import install


def get_openmp():
    """Try to compile/link an example program to check for OpenMP support.
    
    Based on:
    1) http://stackoverflow.com/questions/16549893/programatically-testing-for-openmp-support-from-a-python-setup-script
    2) https://github.com/lpsinger/healpy/blob/6c3aae58b5f3281e260ef7adce17b1ffc68016f0/setup.py
    """
    
    import shutil
    from distutils import sysconfig
    from distutils import ccompiler
    compiler = ccompiler.new_compiler()
    sysconfig.get_config_vars()
    sysconfig.customize_compiler(compiler)
    cc = compiler.compiler
    
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)
    
    fh = open('test.c', 'w')
    fh.write(r"""#include <omp.h>
#include <stdio.h>
int main(void) {
#pragma omp parallel
printf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());
return 0;
}
""")
    fh.close()
    
    ccmd = []
    ccmd.extend( cc )
    ccmd.extend( ['-fopenmp', 'test.c', '-o test'] )
    if os.path.basename(cc[0]).find('gcc') != -1:
        ccmd.append( '-lgomp' )
    elif os.path.basename(cc[0]).find('clang') != -1:
        ccmd.extend( ['-L/opt/local/lib/libomp', '-lomp'] )
    try:
        output = subprocess.check_call(ccmd)
        outCFLAGS = ['-fopenmp',]
        outLIBS = []
        if os.path.basename(cc[0]).find('gcc') != -1:
            outLIBS.append( '-lgomp' )
        elif os.path.basename(cc[0]).find('clang') != -1:
            outLIBS.extend( ['-L/opt/local/lib/libomp', '-lomp'] )
            
    except subprocess.CalledProcessError:
        print("WARNING:  OpenMP does not appear to be supported by %s, disabling" % cc[0])
        outCFLAGS = []
        outLIBS = []
        
    finally:
        os.chdir(curdir)
        shutil.rmtree(tmpdir)
        
    return outCFLAGS, outLIBS


class dummy_install(install):
    """Dummay install method that doesn't let you install."""
    def finalize_options(self, *args, **kwargs):
        raise RuntimeError("This is a dummy package that cannot be installed")


openmpFlags, openmpLibs = get_openmp()


coreExtraFlags = openmpFlags
coreExtraFlags.append('-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION')
coreExtraLibs = openmpLibs


ExtensionModules = [Extension('_helper', ['helper.c',],
                              include_dirs=[numpy.get_include()], libraries=['m'],
                              extra_compile_args=coreExtraFlags, extra_link_args=coreExtraLibs)]


setup(
    cmdclass = {'install': dummy_install}, 
    name = 'dummy_package',
    version = '0.0',
    description = 'This is a dummy package to help build the pulsar extensions',
    ext_modules = ExtensionModules
)
