
ifdef USE_INT
MACRO = -DUSE_INT
endif

#compiler setup
PCC = g++
PCFLAGS = -std=c++14 -O3 -pthread $(MACRO)

COMMON= core/utils.h core/cxxopts.h core/get_time.h core/graph.h core/quick_sort.h

SERIAL= triangle_counting page_rank
PARALLEL= triangle_counting_parallel page_rank_parallel
ALL= $(SERIAL) $(PARALLEL)


all : $(ALL)

% : %.cpp $(COMMON)
	$(PCC) $(PCFLAGS) -o $@ $<

.PHONY : clean

clean :
	rm -f *.o *.obj $(ALL)
