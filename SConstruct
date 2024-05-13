import os, sys, glob, json

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

## Source files (for lib)
lib_src_files = glob.glob(r"src/*.cpp")
lib_src_files += glob.glob(r"src/icgem/*.cpp")
lib_src_files += glob.glob(r"src/free_core_nutation/*.cpp")
lib_src_files += glob.glob(r"src/iau/*.cpp")
lib_src_files += glob.glob(r"src/aod1b/*.cpp")
lib_src_files += glob.glob(r"src/doodson/*.cpp")
lib_src_files += glob.glob(r"src/solid_earth_tides/*.cpp")
lib_src_files += glob.glob(r"src/atmospheric_tides/*.cpp")
lib_src_files += glob.glob(r"src/pole_tide/*.cpp")
lib_src_files += glob.glob(r"src/earth_rotation/*.cpp")
lib_src_files += glob.glob(r"src/eop/*.cpp")
lib_src_files += glob.glob(r"src/gravity/*.cpp")
lib_src_files += glob.glob(r"src/planets/*.cpp")
lib_src_files += glob.glob(r"src/fundarg/*.cpp")
lib_src_files += glob.glob(r"src/relativity/*.cpp")
#lib_src_files += glob.glob(r"src/grid/*.cpp")
#lib_src_files += glob.glob(r"src/trp/*.cpp")

## Headers (for lib)
hdr_src_files = glob.glob(r"src/*.hpp")

## Environments ...
denv = Environment(PREFIX=GetOption(
    'prefix'), CXXFLAGS='-std=c++17 -g -pg -Wall -Wextra -Werror -pedantic -W -Wshadow -Winline -Wdisabled-optimization -DDEBUG -DEIGEN_NO_AUTOMATIC_RESIZING')
penv = Environment(PREFIX=GetOption(
    'prefix'), CXXFLAGS='-std=c++17 -Wall -Wextra -Werror -pedantic -W -Wshadow -Winline -O2 -march=native -DEIGEN_NO_AUTOMATIC_RESIZING')

## Command line arguments ...
debug = ARGUMENTS.get('debug', 0)
test = ARGUMENTS.get('test', 0)
sofa = ARGUMENTS.get('test-vs-sofa', 0)
costg = ARGUMENTS.get('test-vs-costg', 0)
examples = ARGUMENTS.get('make-examples', 0)

## Construct the build enviroment
env = denv.Clone() if int(debug) else penv.Clone()

## What compiler should we be using ?
if GetOption('cxx') is not None: env['CXX'] = GetOption('cxx')

## Set the C++ standard
cxxstd = GetOption('std')
env.Append(CXXFLAGS=' --std=c++{}'.format(cxxstd))

## (shared) library ...
vlib = env.SharedLibrary(source=lib_src_files, target=lib_name, CPPPATH=[
                         'src/'], SHLIBVERSION=lib_version)

## Build ....
#env.Program(source='src/hardisp.cpp',
#    target='bin/hardisp',
#    LIBS=vlib+['geodesy', 'datetime', 'sofa_c'], 
#    LIBPATH='.', 
#    CPPPATH=['src/'])
env.Alias(target='install', source=env.Install(dir=os.path.join(
    GetOption('prefix'), 'include', inc_dir), source=hdr_src_files))
env.Alias(target='install', source=env.InstallVersionedLib(
    dir=os.path.join(GetOption('prefix'), 'lib'), source=vlib))

## Tests ...
if test:
  print('--> Compiling tests ...')
  tenv = env.Clone()
  tenv['CXXFLAGS'] = ' '.join([ x for x in env['CXXFLAGS'].split() if 'inline' not in x])
  cmp_error_fn = 'test/unit_tests/test_compilation_error.json'
  cerror_dct = {}
  if os.path.isfile(cmp_error_fn): os.remove(cmp_error_fn)
  tests_sources = glob.glob(r"test/unit_tests/*.cpp")
  #tests_sources += glob.glob(r"test/time/*.cpp")
  tenv.Append(RPATH=root_dir)
  for tsource in tests_sources:
    ttarget = os.path.join(os.path.dirname(tsource), os.path.basename(tsource).replace('_', '-').replace('.cpp', '.out'))
    if 'mock' in os.path.basename(tsource):
      cerror_dct[os.path.basename(tsource)] = {
                'name': '{:}'.format(os.path.abspath(tsource)),
                'cxx': '{:}'.format(tenv['CXX']),
                'incp' : '{:}'.format(os.path.abspath(os.path.join(tenv['RPATH'], 'src'))),
                'flags': '{:}'.format(' '.join(['-o', tenv['CXXFLAGS']])), 
                'exit': 1}
    else:
      tenv.Program(target=ttarget, source=tsource, CPPPATH='src/', LIBS=vlib+['geodesy', 'datetime', 'cspice'], LIBPATH='.')
    with open(cmp_error_fn, 'w') as fjson: print(json.dumps(cerror_dct, indent = 4), file=fjson)

## SOFA Tests ... 
if sofa:
    print('Compiling tests against SOFA')
    tenv = env.Clone()
    tenv['CXXFLAGS'] = ' '.join([ x for x in env['CXXFLAGS'].split() if 'inline' not in x])
    test_sources = glob.glob(r"test/sofa/*.cpp")
    tenv.Append(RPATH=root_dir)
    for tsource in test_sources:
        ttarget = os.path.join(os.path.dirname(tsource), os.path.basename(tsource).replace('_', '-').replace('.cpp', '.out'))
        tenv.Program(target=ttarget, source=tsource, CPPPATH='src/', LIBS=vlib+['sofa_c', 'geodesy', 'datetime', 'cspice'], LIBPATH='.')

## Compile examples
if examples:
    print('Compiling examples ...')
    tenv = env.Clone()
    tenv['CXXFLAGS'] = ' '.join([ x for x in env['CXXFLAGS'].split() if 'inline' not in x])
    test_sources = glob.glob(r"examples/*.cpp")
    tenv.Append(RPATH=root_dir)
    for tsource in test_sources:
        ttarget = os.path.join(os.path.dirname(tsource), os.path.basename(tsource).replace('_', '-').replace('.cpp', '.out'))
        tenv.Program(target=ttarget, source=tsource, CPPPATH='src/', LIBS=vlib+[ 'geodesy', 'datetime', 'cspice'], LIBPATH='.')

## COSTG Tests ... 
if costg:
    print('Compiling tests against COSTG')
    tenv = env.Clone()
    tenv['CXXFLAGS'] = ' '.join([ x for x in env['CXXFLAGS'].split() if 'inline' not in x])
    test_sources = glob.glob(r"costg-benchmark/bin/*.cpp")
    test_sources = [ fn for fn in test_sources if 'costg_utils.cpp' not in fn ] 
    tenv.Append(RPATH=root_dir)
    for tsource in test_sources:
        ttarget = os.path.join(os.path.dirname(tsource), os.path.basename(tsource).replace('_', '-').replace('.cpp', '.out'))
        tenv.Program(target=ttarget, source=[tsource,'costg-benchmark/bin/costg_utils.cpp'], CPPPATH='src/', LIBS=vlib+[ 'geodesy', 'datetime', 'cspice'], LIBPATH='.')
