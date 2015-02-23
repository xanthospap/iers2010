noinst_LIBRARIES                  = libiers10-alternatives.a
libiers10_alternatives_a_CXXFLAGS = -std=c++11 -g -Wall -Wextra -Werror -pedantic -W -Wshadow -Winline -Wdisabled-optimization
dist_libiers10_alternatives_a_SOURCES = gpt.cpp \
	gmf.cpp  \
	gpt.cpp \
	shclass.cpp \
	shclass.hpp
noinst_PROGRAMS = alternatives
alternatives_SOURCES = main.cpp
alternatives_CXXFLAGS = -std=c++11 -g -Wall -Wextra -Werror -pedantic -W -Wshadow -Winline -Wdisabled-optimization
alternatives_LDADD = libiers10-alternatives.a