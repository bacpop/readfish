from codecs import open
from os import path
import os
import re
import sys
import io
import platform
import subprocess

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion

PKG_NAME = "readfish"
MOD_NAME = "ru"

DESCRIPTION = """# <img src="https://raw.githubusercontent.com/LooseLab/ru/rename_cli/examples/images/readfish_logo.jpg">

Installation
---

This toolkit currently requires MinKNOW (minknow-core v4.0.4) to be installed and 
[`read_until_api`](https://github.com/nanoporetech/read_until_api) to be installed
separately. We recommend installing in a virtual environment as so:

```bash
# Make a virtual env
python3 -m venv readfish
source ./readfish/bin/activate
pip install git+https://github.com/nanoporetech/read_until_api
pip install readfish
```

Usage
---
```bash
# check install
$ readfish
usage: readfish [-h] [--version]
                {targets,align,centrifuge,unblock-all,validate,summary} ...

positional arguments:
  {targets,align,centrifuge,unblock-all,validate,summary}
                        Sub-commands
    targets             Run targeted sequencing
    align               ReadFish and Run Until, using minimap2
    centrifuge          ReadFish and Run Until, using centrifuge
    unblock-all         Unblock all reads
    validate            ReadFish TOML Validator
    summary             Summary stats from FASTQ files

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit

See '<command> --help' to read about a specific sub-command.

# example run command - change arguments as necessary:
$ readfish targets --experiment-name "Test run" --device MN17073 --toml example.toml --log-file RU_log.log
```
"""

def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")
    ) as fp:
        return fp.read()

def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable,
                      '-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON']

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j1']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)

here = path.abspath(path.dirname(__file__))

with open(path.join(here, "README.md")) as fh, open(
    path.join(here, "requirements.txt")
) as req:
    install_requires = [pkg.strip() for pkg in req]

__version__ = ""
exec(open("{}/_version.py".format(MOD_NAME)).read())

setup(
    name=PKG_NAME,
    version=__version__,
    author="Alexander Payne",
    author_email="alexander.payne@nottingham.ac.uk",
    description="Adaptive sampling toolkit for MinION",
    long_description=DESCRIPTION,
    long_description_content_type="text/markdown",
    url="https://github.com/LooseLab",
    packages=find_packages(exclude=["*.test", "*.test.*", "test.*", "test"]),
    entry_points={
        "console_scripts": [
            "ru_validate={}.validate:main".format(MOD_NAME),
            "ru_generators={}.ru_gen:main".format(MOD_NAME),
            "ru_summarise_fq={}.summarise_fq:main".format(MOD_NAME),
            "ru_iteralign={}.iteralign:main".format(MOD_NAME),
            "ru_iteralign_centrifuge={}.iteralign_centrifuge:main".format(MOD_NAME),
            "ru_unblock_all={}.unblock_all:main".format(MOD_NAME),
            "readfish={}.cli:main".format(MOD_NAME),
            "ru_generate_graph={}.graphalign.generate_graph:main".format(MOD_NAME),
        ],
    },
    ext_modules=[CMakeExtension('query_cpp')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    install_requires=install_requires,
    include_package_data=True,
    python_requires=">=3.5",
)
