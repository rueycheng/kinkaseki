bin_PROGRAMS = tokenize glue evaluate_segmentation ubt descent

tokenize_SOURCES = tokenize.cpp 
glue_SOURCES = glue.cpp compat.h 
evaluate_segmentation_SOURCES = evaluate_segmentation.cpp 
ubt_SOURCES = ubt.cpp compat.h
descent_SOURCES = descent.cpp

LDADD = -lboost_program_options -lboost_filesystem -lboost_regex -lboost_system
AM_CXXFLAGS = -Wall -std=c++0x -I@topdir@
