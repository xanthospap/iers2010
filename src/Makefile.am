lib_LTLIBRARIES = libiers2010.la

## Only use the following flags for debuging purposes
libiers2010_la_CXXFLAGS = \
	-std=c++14 \
	-g \
	-pg \
	-Wall \
	-Wextra \
	-Werror \
	-pedantic \
	-W \
	-Wshadow \
	-Winline \
	-Wdisabled-optimization \
	-DDEBUG

dist_include_HEADERS = \
	iers2010.hpp

dist_libiers2010_la_SOURCES = \
	fundarg.cpp \
	pmsdnut2.cpp \
	utlibr.cpp
