CPP=g++
CPPFLAGS=-g -Wall -pthread -std=c++17
LIBS=-lz -lboost_program_options -ltbb
ifneq ($(DEBUG),)
	CPPFLAGS+=-O0 -DDBGPRINT
else
	CPPFLAGS+=-O2
endif

HEADERS = $(wildcard *.h)
SRCS = $(wildcard *.cpp) ksw2_extz2_sse.c
OBJS = $(SRCS:.cpp=.o)

.PHONY: all
all: motif-search

motif-search: $(OBJS)
	$(CPP) $(CPPFLAGS) -o $@ $(OBJS) ${LIBS}

sample: motif-search
	./motif-search -r `ls sample_data/reads/EINS3-guppy2/*.fastq | xargs echo -n` -m sample_data/motifs/EINS.fa -l 5

sample-query: motif-search
	./motif-search -r sample_data/multiple_queries/merge.fa -m sample_data/multiple_queries/256Output/all5P-1.fa -q sample_data/multiple_queries/256Output/all3P.fa

clean:
	rm motif-search *.o
