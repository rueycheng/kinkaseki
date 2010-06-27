bin_PROGRAMS = tokenize collect compile kmeans blm
tokenize_SOURCES = tokenize.cpp compat.h magicbox.h porter.c porter.h
collect_SOURCES = collect.cpp compat.h magicbox.h
compile_SOURCES = compile.cpp compat.h magicbox.h
kmeans_SOURCES = kmeans.cpp compat.h magicbox.h
blm_SOURCES = blm.cpp compat.h magicbox.h
#--------------------------------------------------
# bin_PROGRAMS = gii_parse gii_tokenize gii_index gii_model1 gii_model2 gii_model3 sim merge lda slda
# gii_parse_SOURCES = gii_parse.cpp magicbox.h expatpp.h expatpp_callback.h gii_parse.h
# gii_tokenize_SOURCES = gii_tokenize.cpp porter.c porter.h magicbox.h gii_tokenize.h
# gii_index_SOURCES = gii_index.cpp magicbox.h
# gii_model1_SOURCES = gii_model1.cpp magicbox.h porter.h porter.c
# gii_model2_SOURCES = gii_model2.cpp magicbox.h porter.h porter.c
# gii_model3_SOURCES = gii_model3.cpp magicbox.h porter.h porter.c
# sim_SOURCES = sim.cpp magicbox.h
# merge_SOURCES = merge.cpp magicbox.h
# lda_SOURCES = lda.cpp magicbox.h
# slda_SOURCES = slda.cpp magicbox.h
#-------------------------------------------------- 

#--------------------------------------------------
# bin_PROGRAMS += prep
# prep_SOURCES = prep.cpp magicbox.h
#-------------------------------------------------- 

AM_LDFLAGS = -lboost_program_options-mt -lboost_filesystem-mt
AM_CXXFLAGS = -Wall

if USE_CXX0X
    AM_CXXFLAGS += -std=c++0x
endif

INCLUDES = -I@topdir@
