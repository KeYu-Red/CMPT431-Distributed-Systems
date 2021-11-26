#compiler setup
CXX = g++
CXXFLAGS = -std=c++14 -march=native -pthread -O3 
LDFLAGS = -Llib -lalloc431 -lpthread

THROUGHPUT=one_lock_queue_throughput two_lock_queue_throughput non_blocking_queue_throughput 
CORRECTNESS=one_lock_queue_correctness two_lock_queue_correctness non_blocking_queue_correctness 
ALL= $(THROUGHPUT) $(CORRECTNESS)
ONE_LOCK_QUEUE = queues/one_lock_queue.h
TWO_LOCK_QUEUE = queues/two_lock_queue.h
NON_BLOCKING_QUEUE = queues/non_blocking_queue.h
COMMON= common/utils.h common/cxxopts.h common/get_time.h common/quick_sort.h common/parallel.h common/allocator.h

all : $(ALL)

# ------------ THROUGHPUT --------------
one_lock_queue_throughput : driver_throughput.cpp $(COMMON) $(ONE_LOCK_QUEUE)
	$(CXX) $(CXXFLAGS) -DONE_LOCK_QUEUE $< $(LDFLAGS) -o $@

two_lock_queue_throughput : driver_throughput.cpp $(COMMON) $(TWO_LOCK_QUEUE)
	$(CXX) $(CXXFLAGS) -DTWO_LOCK_QUEUE $< $(LDFLAGS) -o $@

non_blocking_queue_throughput : driver_throughput.cpp $(COMMON) $(NON_BLOCKING_QUEUE)
	$(CXX) $(CXXFLAGS) -DNON_BLOCKING_QUEUE $< $(LDFLAGS) -o $@

# ------------ CORRECTNESS --------------
one_lock_queue_correctness : driver_correctness.cpp $(COMMON) $(ONE_LOCK_QUEUE)
	$(CXX) $(CXXFLAGS) -DONE_LOCK_QUEUE $< $(LDFLAGS) -o $@

two_lock_queue_correctness : driver_correctness.cpp $(COMMON) $(TWO_LOCK_QUEUE)
	$(CXX) $(CXXFLAGS) -DTWO_LOCK_QUEUE $< $(LDFLAGS) -o $@

non_blocking_queue_correctness :driver_correctness.cpp $(COMMON) $(NON_BLOCKING_QUEUE)
	$(CXX) $(CXXFLAGS) -DNON_BLOCKING_QUEUE $< $(LDFLAGS) -o $@

generate_test_input : generate_test_input.cpp $(COMMON)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(COMMON):

.PHONY : clean

clean :
	rm -f *.o *.obj $(ALL)

cleansrc :
	rm -f *.o $(ALL) $(OTHERS)

