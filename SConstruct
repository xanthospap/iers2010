import os, sys, glob

## Prefic for install(ed) files
## prefix="/usr/local"
## Library version
lib_version="0.1.0"
## Library name
lib_name="iers2010"
## Include dir (following prefix) if any
inc_dir="iers2010"
## the rootdir of the project
root_dir=os.path.abspath(os.getcwd())

## get number of CPUs and use for parallel builds
num_cpu = int(os.environ.get('NUM_CPU', 2))
SetOption('num_jobs', num_cpu)
print("running with -j %s" % GetOption('num_jobs'))

## user can specify the --prefix variables (expanded to $PREFIX)
AddOption('--prefix',
          dest='prefix',
          type='string',
          nargs=1,
          action='store',
          metavar='DIR',
          help='installation prefix',
          default='/usr/local')
AddOption('--cxx',
          dest='cxx',
          type='string',
          nargs=1,
          action='store',
          metavar='CXX',
          help='C++ Compiler',
          default=None)
AddOption('--std',
          dest='std',
          type='string',
          nargs=1,
          action='store',
          metavar='STD',
          help='C++ Standard [11/14/17/20]',
          default='17')
AddOption('--make-check',
          dest='check',
          action='store_true',
          help='Trigger building of test programs',
          default=False)
AddOption('--make-sofa-check',
          dest='sofa_check',
          action='store_true',
          help='Trigger building of test programs against the SOFA IAU library. The library must be available in the system (header/library)',
          default=False)

## Source files (for lib)
lib_src_files = glob.glob(r"src/*.cpp")
lib_src_files += glob.glob(r"src/iau/*.cpp")
lib_src_files += glob.glob(r"src/ch5/*.cpp")
lib_src_files += glob.glob(r"src/ch7/*.cpp")
lib_src_files += glob.glob(r"src/ch8/*.cpp")
lib_src_files += glob.glob(r"src/ch9/*.cpp")
lib_src_files += glob.glob(r"src/hardisp/*.cpp")
lib_src_files += glob.glob(r"src/dehanttideinel/*.cpp")
lib_src_files += glob.glob(r"src/extra/atmosphere/*.cpp")
lib_src_files += glob.glob(r"src/interpf/*.cpp")

## Headers (for lib)
hdr_src_files = glob.glob(r"src/*.hpp")

## Environments ...
denv = Environment(PREFIX=GetOption(
    'prefix'), CXXFLAGS='-std=c++17 -g -pg -Wall -Wextra -Werror -pedantic -W -Wshadow -Winline -Wdisabled-optimization -DDEBUG')
penv = Environment(PREFIX=GetOption(
    'prefix'), CXXFLAGS='-std=c++17 -Wall -Wextra -Werror -pedantic -W -Wshadow -Winline -O2 -march=native')

## Command line arguments ...
debug = ARGUMENTS.get('debug', 0)
eigen = ARGUMENTS.get('eigen', 0)

## Construct the build enviroment
env = denv.Clone() if int(debug) else penv.Clone()

## What compiler should we be using ?
if GetOption('cxx') is not None: env['CXX'] = GetOption('cxx')

## Use eigen ?
if eigen: env.Append(CXXFLAGS=' -DUSE_EIGEN')

## Set the C++ standard
cxxstd = GetOption('std')
env.Append(CXXFLAGS=' --std=c++{}'.format(cxxstd))

## (shared) library ...
vlib = env.SharedLibrary(source=lib_src_files, target=lib_name, CPPPATH=[
                         'src/'], SHLIBVERSION=lib_version)

## Build ....
env.Program(source='src/hardisp.cpp', target='bin/hardisp',
            LIBS=vlib+['geodesy', 'datetime','matvec'], LIBPATH='.', CPPPATH=['src/'])
env.Alias(target='install', source=env.Install(dir=os.path.join(
    GetOption('prefix'), 'include', inc_dir), source=hdr_src_files))
env.Alias(target='install', source=env.InstallVersionedLib(
    dir=os.path.join(GetOption('prefix'), 'lib'), source=vlib))

## Tests ...
if GetOption('check') is not None and GetOption('check'):
  tests_sources = glob.glob(r"test/*.cpp")
  env.Append(RPATH=root_dir)
  for tsource in tests_sources:
    pth = os.path.dirname(tsource)
    bsn = os.path.basename(tsource)
    ttarget = os.path.join(pth, bsn.replace(
        '_', '-').replace('.cpp', '.out'))
    env.Program(target=ttarget, source=tsource, CPPPATH='src/',
                LIBS=vlib+['geodesy', 'datetime','matvec'], LIBPATH='.')

## Build tests against SOFA lib
if GetOption('sofa_check') is not None and GetOption('sofa_check'):
  tests_sources = glob.glob(r"test/sofa_unit_tests/sofa*.cpp")
  env.Append(RPATH=root_dir)
  for tsource in tests_sources:
    pth = os.path.dirname(tsource)
    bsn = os.path.basename(tsource)
    ttarget = os.path.join(pth, bsn.replace(
        '_', '-').replace('.cpp', '.out'))
    env.Program(target=ttarget, source=tsource, CPPPATH='src/',
                LIBS=vlib+['geodesy', 'datetime','matvec', 'sofa_c'], LIBPATH='.')
