import os, sys, glob

## Prefic for install(ed) files
PREFIX="/usr/local"
if not os.path.isdir(PREFIX):
    print('[ERROR] Cannot find \'prefix\' directory, aka {:}; aborting'.format(PREFIX), file=sys.stderr)
    sys.exit(1)

## Library version
LIBV="0.1.0"
## Library name
LIBNAME="iers2010"
## Include dir (following prefix) if any
INCDIR="iers2010"

## get number of CPUs and use for parallel builds
num_cpu = int(os.environ.get('NUM_CPU', 2))
SetOption('num_jobs', num_cpu)
print("running with -j %s" % GetOption('num_jobs'))

## Source files (for lib)
lib_src_files = glob.glob(r"src/*.cpp")

## Headers (for lib)
hdr_src_files = glob.glob(r"src/*.hpp")

## Environments ...
denv = Environment(CXXFLAGS='-std=c++17 -g -pg -Wall -Wextra -Werror -pedantic -W -Wshadow -Winline -Wdisabled-optimization -DDEBUG')
penv = Environment(CXXFLAGS='-std=c++17 -Wall -Wextra -Werror -pedantic -W -Wshadow -Winline -O2 -march=native')

## Command line arguments ...
debug = ARGUMENTS.get('debug', 0)

## Construct the build enviroment
env = denv.Clone() if int(debug) else penv.Clone()

## (shared) library ...
vlib = env.SharedLibrary(source=lib_src_files, target=LIBNAME, CPPPATH=['.'], SHLIBVERSION=LIBV)

## Build ....
env.Program(source='src/hardisp.cpp', target='bin/hardisp',
            LIBS=vlib+['ggeodesy', 'ggdatetime'], LIBPATH='.', CPPPATH=['src/'])
env.Alias(target='install', source=env.Install(dir=os.path.join(PREFIX, 'lib'), source=vlib))
env.Alias(target='install', source=env.Install(dir=os.path.join(PREFIX, 'include', INCDIR), source=hdr_src_files))

## Tests ...
tests_sources = glob.glob(r"test/*.cpp")
env.Append(CXXFLAGS='-Wl,-rpath,\'\'')
for tsource in tests_sources:
  ttarget = tsource.replace('_', '-').replace('.cpp', '.out')
  env.Program(target=ttarget, source=tsource, CPPPATH='src/', LIBS=vlib+['ggeodesy', 'ggdatetime'], LIBPATH='.')
