**********************************************
(c)2016: Oxford Nanopore Technologies, Limited
**********************************************


Nanonet
=======

Nanonet provides recurrent neural network basecalling for Oxford Nanopore MinION
data. It represents the first generation of such a basecaller from Oxford Nanopore
Technologies, and is provided as a technology demonstrator. Nanonet is provided
unsupported by Oxford Nanopore Technologies, see LICENSE.md for more information.

For training networks, Nanonet leverages currennt to run recurrent neural
networks. Currennt is generally run with GPUs to aid performance but can be run
in a CPU only environment. The basecaller does not require currennt, and is
written in pure python with minimal requirements.


Installation
------------

*Nanonet contains implementations of both 1D and 2D basecalling, with OpenCL
versions of each of these. By default, using the instructions in this section
only the canonical 1D basecalling library will be compiled; OpenCL acceleration
of 1D basecalling and any 2D basecalling support will not be configured. See later
sections of this documented for setting up these components.*

**Requirements**

Nanonet is a python-based commandline suite consisting of several programs for
performing basecalling with recurrent neural networks. It has been developed
using 64-bit Python 2.7 available from python.org.

In addition to Python the only required dependencies are h5py and numpy. These can
be downloaded and installed/compiled automatically, though installing them from your
system's package repository is generally preferable in the first instance. For Windows
Christophe Golke maintains a repository of compiled python wheels at:

    http://www.lfd.uci.edu/~gohlke/pythonlibs/

For OSX homebrew is easiest:

    brew tap homebrew/science
    brew install hdf5

**Easy Windows Installation**

Python wheels for nanonet are available from the github page under the releases
tab:

    https://github.com/nanoporetech/nanonet-dev/releases

to install these run the following from a command prompt:

    pip install nanonet-2.0.0-cp27-cp27m-win_amd64.whl

The wheel currently enables CPU implementataions of both 1D and 2D calling,
and additionally the OpenCL implementation of 1D calling. To enable the latter
please refer to the OpenCL section below. Wheels for other platforms may be made
available in the future.

**Setup**

The basecalling component of nanonet should install quite trivially using the
standard python mechanism on most platforms:

    python setup.py install

The basecaller contains small amounts of C code for performing event detection.
A C compiler is therefore required. Under Linux and OSX, one is likely
installed already, for Windows one can download and install the Microsoft Visual
C++ Compiler for Python 2.7 from:

    https://www.microsoft.com/en-gb/download/details.aspx?id=44266

Most of these documentation assume you are using this compiler on Windows. The
section **Compiling with MinGW-w64** explains how to use GCC based compilation
on Windows.

**Optional Watcher Component**

Nanonet contains an optional component to watch a filesystem for read files as
they are produced by MinKnow. This feature is not installed by default. To
install it run

    pip install -e  .[watcher]

from the source directory, or simply

    pip install watchdog

for any location. This will allow use of the `--watch` option of the basecaller.


Peforming basecalling
---------------------

Nanonet provides a single program for basecalling Oxford Nanopore Technolgies'
reads from .fast5 files. The output is to stdout as fasta formatted sequences.

If you followed the instructions above the program `nanonetcall` should be on
your path. It requires a single argument:

    nanonetcall {input_folder} > {output.fasta}

To test your installation several .fast5 files are provided in the
`sample_data` folder of the source distribution, as a concrete example the
command:

    nanonetcall --jobs 5 --chemistry r9 sample_data > basecalls.fa

will produced output along the following lines:

    Basecalled 5 reads (25747 bases, 49340 events) in 36.6668331623s (wall time)
    Profiling
    ---------
    Feature generation: 1.08410215378
    Load network: 0.0128238201141
    Run network: 15.6343309879 (1.5990450771 kb/s, 3.13412835112 kev/s)
    Decoding: 19.8838129044 (1.25730412574 kb/s, 2.46431608644 kev/s)
    
**Filesystem watching**

Nanonet has the ability to watch a filesystem as reads are produced. This
behaviour is enabled with the `--watch` option:

    nanonetcall input_folder --watch 600 > basecalls.fa

where the option value is a timeout in seconds, when no new reads are seen for
this time nanonetcall will exit.

**Input files**

nanonetcall operates from single-read .fast5 files as output by MinKnow. These
should contain raw data; the event detection step and segmentation into template
and complement sections will be performed by nanonet.

**Using multiple CPUs**

By default nanonetcall will use a maximum of one CPU. This can be altered
through use of the `--jobs` option. In using this option be aware that
higher numbers will lead to increased memory usage.


2D Basecalling
--------------

The 2D basecalling library has a dependency on the boost C++ library. For this
reason it is not compiled by default. If you know that you have a working boost
installation you can simply use the following to enable 2D basecalling:

    python setup.py install with2d
    
Boost can be installed on most Linux systems via the system package manager. For
example on Ubuntu it should be sufficient to install the boost-python package:

    sudo apt-get install libboost-python1.54-dev libboost-python1.54.0

You may elect to use a different version of boost if you wish. On OSX boost can
again use homebrew:

    brew install boost --with-python
    brew install boost-python
    
On Windows the simplest method to install boost is to obtain a precompiled version
from sourceforge:

    https://sourceforge.net/projects/boost/files/boost-binaries/

nanonet has been tested with version 1.55.0 on windows. It is important to use the
version compiled with the same version of the Microsoft Visual C compiler you are
using. If you followed the instructions above to install the Microsoft Visual
C++ Compiler for Python 2.7 you should download the package labelled `msvc-9.0-64.exe`,
i.e. the file available here:

    https://sourceforge.net/projects/boost/files/boost-binaries/1.55.0/boost_1_55_0-msvc-9.0-64.exe/download
    
The above installer will by default install boost to `c:\local\boost_1_55_0`. If
you change this path you will also need to edit the `setup.py` file in nanonet. To run
the 2D basecaller you will need also to add the library location to your `Path` environment
variable. In Windows Powershell this can be done with:

    $env:Path += ";c:\local\boost_1_55_0\lib64-msvc-9.0"

Once you have installed boost on your OS, the 2D basecalling components can be
compiled and set up with:

    python setup.py install with2d
    
Performing 2D basecalling currently requires use of a distinct program from the pure 1D
basecaller. The interface of this program is much the same as `nanonetcall`, for example
a basic use would simply require:

    nanonet2d --chemistry r9 sample_data calls
    
The second option here specifies a prefix for output fasta files; three files will be
created: one each for template, complement and 2D basecalls.


OpenCL Support
--------------

Nanonet contains OpenCL accelerated versions of both 1D and 2D basecalling. Currently
the implementations of these use different mechanisms for creating the OpenCL kernels:
1D basecalling acceleration is driven through `pyopencl` whilst 2D acceleration directly
interfaces with OpenCL libraries.

Configuring a working OpenCL environment depends heavily on your OS and device you wish
to target. If you are unfamiliar with how to do this you are recommended to start with
the `pyopencl` documentation:

    https://wiki.tiker.net/PyOpenCL

**1D Acceleration**

Once you have the `pyopencl` examples working you should be able to run the accelerated
version of 1D basecalling:

    nanonetcall --chemistry r9 --platforms <VENDOR:DEVICE:1> --exc_opencl sample_data

where `VENDOR` and `DEVICE` will be machine dependent. To see available devices you can
examine the output of:

    nanonetcall --list_platforms
    
The `--exc_opencl` option above instructs nanonet to use only the OpenCL device(s) listed
on the command line. Without this option, CPU resources will also be used. You may wish to
experiment with this option and the `--jobs` option to achieve optimal throughput.

**2D Acceleration**

To enable OpenCL acceration on 2D basecalling it is necessary to compile additional
components during the setup of nanonet by using the `opencl2d` argument:

    python setup.py install opencl2d

With a working OpenCL runtime and development environment 2D basecalling can be accelerated
by simply adding and option on the commandline:

    nanonet2d --opencl_2d --chemistry r9 sample_data calls
    
The program will automatically choose an OpenCL device to use, giving preference
to GPU devices over CPU ones. It is not currently possible to use OpenCL acceleration for
the 1D basecalling necessary for performing a 2D basecall.

Compiling with MinGW-w64
------------------------

The Microsoft Visual C++ Compiler for Python 2.7 compiler is fairly old and can
produce less than efficient code compared to more recent compilers. The
following describes how one may instead compile nanonet using MinGW. This process
is not recommended for those unfamiliar with compilation toolchains. We use
specific, but not the latest, versions of MinGW and Boost in what follows. Be
aware that these instructions may not work with later versions of these tools.

First download the following build of MinGW-w64:

    https://sourceforge.net/projects/mingwbuilds/files/host-windows/releases/4.8.1/64-bit/threads-posix/seh/x64-4.8.1-release-posix-seh-rev5.7z/download

together with Boost 1.55.0:

    https://sourceforge.net/projects/boost/files/boost/1.55.0/boost_1_55_0.7z/download
    
The files are packaged with 7zip, which you should install if you haven't
previously. Extract both of the packages to:

    C:\local\

Open Windows Powershell and run the following to build the boost-python library:

    $env:Path += ';C:\local\mingw64\bin'
    cd c:\local\boost_1_55_0\
    .\bootstrap.bat gcc
    .\b2.exe toolset=gcc address-model=64 link=shared define=MS_WIN64 --with-python stage
    
The `distutils` module of Python 2.7 supports using MinGW as a compiler,
although unfortunately it is likely to cause .dlls to be incorrectly
linked when using MinGW-w64. For this reason and to streamline the build
process, nanonet patches the class `distutils.cygwinccompiler.Mingw32CCompiler`
with its own version and manipulates some `distutils` internals in order
to force use of this class when a users elects to use MinGW-w64. All this
means is that to compile and setup an in-place (development) install of
nanonet run:

    python setup.py develop --user with2d mingw

If you have and existing MinGW-w64 setup you can use the environment
variables `BOOST_ROOT`, `BOOST_LIB`, and `BOOST_PYTHON` to specify
respectively the location of you boost install, relative location of
the boost libraries, and the name of your boost-python library. For
example the defaults are:

    BOOST_ROOT  = "c:\local\boost_1_55_0\"
    BOOST_LIB   = "stage\lib"
    BOOST_PYTON = "boost_python-mgw48-mt-1_55"


Benchmarks
----------

The results below are indicative and provided for reference only. The
performance of nanonet should not be assumed to be equal to that of
a production implementation.

**AWS g2.2xlarge Ubuntu 16.04 Xeon E5-2670 @ 2.60GHz**
- numpy is linked to a multithreaded blas, it appears
  best to set the number of threads to 1.
- 1D performance scales nicely with number of CPU processes
- 1D GPU performance is rough 4.5X for single concurrent kernel
  or approaching 6.5X for ten concurrent kernels.
- 2D GPU performance is 7X.

*1D basecalling*

    $ nanonetcall more_sample_data/ > /dev/null
    Basecalled 20 reads (30160 bases, 62600 events) in 39.9055719376s (wall time)
    Run network:  23.59 ( 1.278 kb/s,  2.653 kev/s)
    Decoding:     14.53 ( 2.076 kb/s,  4.308 kev/s)
    $
    $ OMP_NUM_THREADS=1 nanonetcall more_sample_data/ > /dev/null
    Basecalled 20 reads (30160 bases, 62600 events) in 30.2339270115s (wall time)
    Run network:  13.78 ( 2.189 kb/s,  4.544 kev/s)
    Decoding:     14.74 ( 2.046 kb/s,  4.247 kev/s)
    $
    $ OMP_NUM_THREADS=1 nanonetcall more_sample_data/ --jobs 4 > /dev/null
    Basecalled 20 reads (30160 bases, 62600 events) in 7.72924900055s (wall time)
    Run network:  13.84 ( 2.180 kb/s,  4.525 kev/s)
    Decoding:     15.00 ( 2.011 kb/s,  4.174 kev/s)
    $
    $ OMP_NUM_THREADS=1 nanonetcall more_sample_data/ --platforms nvidia:0:1 --exc_opencl > /dev/null
    Basecalled 20 reads (30160 bases, 62600 events) in 8.73325586319s (wall time)
    Run network:   3.83 ( 7.868 kb/s, 16.331 kev/s)
    Decoding:      1.54 (19.623 kb/s, 40.730 kev/s)

*2D basecalling*

    $ OMP_NUM_THREADS=1 nanonet2d more_sample_data calls
    Processed 20 reads in 214.572195053s (wall time)
    Template Run network:  14.01 ( 2.153 kb/s,  4.469 kev/s)
    Template Decoding:     15.76 ( 1.914 kb/s,  3.973 kev/s)
    2D calling:           156.85 ( 0.188 kb/s)
    $
    $ OMP_NUM_THREADS=1 nanonet2d more_sample_data calls --opencl_2d
    Processed 20 reads in 78.750191927s (wall time)
    Template Run network:  13.92 ( 2.167 kb/s,  4.498 kev/s)
    Template Decoding:     15.50 ( 1.946 kb/s,  4.038 kev/s)
    2D calling:            21.72 ( 1.355 kb/s)

**Intel NUC5i7RYH Ubuntu 14.04 i7-5557U CPU @ 3.10GHz**
- Again CPU performance better with single threaded blas
- 1D CPU marginally better than AWS machine, though relative
  speed of network and coding altered.
- Two physical cores: performance with 4 processes not
  significantly faster than 2.
- 1D Iris 6100 roughly 3.5X over CPU core. Concurrent kernel
  execution is not handled well.
- 2D GPU 5X over CPU.

*1D basecalling*

    $ nanonetcall more_sample_data/ > /dev/null
    Basecalled 20 reads (30160 bases, 62600 events) in 34.7250430584s (wall time)
    Run network:  16.35 ( 1.845 kb/s,  3.829 kev/s)
    Decoding:     16.26 ( 1.855 kb/s,  3.851 kev/s)
    $
    $ OMP_NUM_THREADS=1 nanonetcall more_sample_data/ > /dev/null
    Basecalled 20 reads (30160 bases, 62600 events) in 28.1761329174s (wall time)
    Run network:  10.19 ( 2.960 kb/s,  6.143 kev/s)
    Decoding:     16.13 ( 1.870 kb/s,  3.881 kev/s)
    $
    $ OMP_NUM_THREADS=1 nanonetcall more_sample_data/ --jobs 2 >/dev/null
    Basecalled 20 reads (30160 bases, 62600 events) in 14.3247690201s (wall time)
    Run network:  10.38 ( 2.907 kb/s,  6.033 kev/s)
    Decoding:     16.11 ( 1.872 kb/s,  3.885 kev/s)
    $
    $ OMP_NUM_THREADS=1 nanonetcall more_sample_data/ --jobs 4 >/dev/null
    Basecalled 20 reads (30160 bases, 62600 events) in 12.9280588627s (wall time)
    Run network:  18.68 ( 1.615 kb/s,  3.352 kev/s)
    Decoding:     29.16 ( 1.034 kb/s,  2.147 kev/s)
    $
    $ OMP_NUM_THREADS=1 nanonetcall more_sample_data/ --platforms intel:0:1 --exc_opencl > /dev/null
    Basecalled 20 reads (30160 bases, 62600 events) in 10.3367981911s (wall time)
    Run network:   5.54 ( 5.444 kb/s, 11.299 kev/s)
    Decoding:      2.88 (10.474 kb/s, 21.739 kev/s)

*2D basecalling*

    $ OMP_NUM_THREADS=1 nanonet2d more_sample_data calls
    Processed 20 reads in 180.179665089s (wall time)
    Template Run network:  10.15 ( 2.970 kb/s,  6.165 kev/s)
    Template Decoding:     16.59 ( 1.818 kb/s,  3.773 kev/s)
    2D calling:           127.52 ( 0.231 kb/s)
    $
    $ OMP_NUM_THREADS=1 nanonet2d more_sample_data calls --opencl_2d > /dev/null
    Processed 20 reads in 80.0777461529s (wall time)
    Template Run network:  10.15 ( 2.973 kb/s,  6.170 kev/s)
    Template Decoding:     16.73 ( 1.803 kb/s,  3.742 kev/s)
    2D calling:            26.19 ( 1.124 kb/s)


**E5-2650 v4 @ 2.20GHz (x2, 24 cores total) GTX 1080 Ubuntu 14.04**
- Timings likely a little off as GPU used by other processes
- Utilise very little (<200Mb) of GPU memory

*1D basecalling*

    $ nanonetcall more_sample_data/ --platforms nvidia:0:20 --exc_opencl > /dev/null
    Basecalled 20 reads (30160 bases, 62600 events) in 3.88349485397s (wall time)
    Run network:   1.32 (22.829 kb/s, 47.384 kev/s)
    Decoding:      0.91 (33.268 kb/s, 69.051 kev/s)

*2D basecalling*

    $ OMP_NUM_THREADS=1 nanonet2d more_sample_data calls --opencl_2d
    Processed 20 reads in 56.9682049751s (wall time)
    Template Run network:  10.10 ( 2.986 kb/s,  6.199 kev/s)
    Template Decoding:     11.60 ( 2.601 kb/s,  5.398 kev/s)
    2D calling:            13.30 ( 2.214 kb/s)
    $
    $ OMP_NUM_THREADS=1 nanonet2d more_sample_data calls --opencl_2d --jobs 4
    Processed 20 reads in 15.9661970139s (wall time)
    Template Run network:  10.19 ( 2.959 kb/s,  6.141 kev/s)
    Template Decoding:     11.74 ( 2.570 kb/s,  5.334 kev/s)
    2D calling:            15.96 ( 1.844 kb/s)



Training a network
------------------

The package provides also an interface to currennt for training networks from
.fast5 files. The type of neural networks implemented by currennt require
labelled input data. **Nanonet does not provide a method for labelling data**.
It does however provide a hook to create labelled data. To install currennt
please refer to the full installation instructions in TRAIN_INSTALL.md.

To run the trainer we specify training data, and validation data. The former
will be used to train the network whilst the latter is used to check that the
network does not become overtrained to the training data.

    nanonettrain --train <training_data> --val <validation_data> \
        --output <output_model_prefix> --model <input_model_spec>

An example input model can be found in `nanonet/data/default_model.tmpl`.
The only other consideration is that the size of the first ("input") layer of
the network must correspond to the feature vectors created by
`nanonet.features.events_to_features`. The nanonettrain program will try to
enforce these considerations. In contructing models one should assign the
input layer a size of `<n_features>` and the final two layers as `<n_states>`,
as in the example.

Training is an intensive process, even on a GPU expect it to take hours not
minutes. It is not recommended to attempt training models without GPU support.


Trouble Shooting
----------------

If you performed a user install (`setup.py install --user`) the `nanonetcall`
program may not be on your path. This is because the location into which
setuptools installs programs for users is not often a default item in the
user's path. On OSX the location is typically:

    /Users/<username>/Library/Python/2.7/bin/

whilst on Ubuntu it is:

    ~/.local/bin

If `nanonet` complains that it cannot locate the `currennt` executable you will
need to set the `CURRENNT` environment variable to the location of the
executable.


