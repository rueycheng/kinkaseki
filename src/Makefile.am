bin_PROGRAMS = tokenize glue evaluate_segmentation

tokenize_SOURCES = tokenize.cpp 
glue_SOURCES = glue.cpp compat.h 
evaluate_segmentation_SOURCES = evaluate_segmentation.cpp 

#--------------------------------------------------
# blm_SOURCES = blm.cpp compat.h magicbox.h
# collect_SOURCES = collect.cpp compat.h magicbox.h
# compile_SOURCES = compile.cpp compat.h magicbox.h
# kmeans_SOURCES = kmeans.cpp compat.h magicbox.h
#-------------------------------------------------- 

AM_LDFLAGS = -lboost_program_options-mt -lboost_filesystem-mt
AM_CXXFLAGS = -Wall -std=c++0x
INCLUDES = -I@topdir@
