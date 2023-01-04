CFLAGS=-g -Wall -O3
CXXFLAGS=$(CFLAGS) -std=c++11
LIBS=-lz
PROG=DBGPS-greedy-path

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.PHONY:all clean

all:$(PROG)

DBGPS-greedy-path:DBGPS-greedy-path.c khashl.h ketopt.h kseq.h kthread.h
	$(CC) $(CFLAGS) -o $@ DBGPS-greedy-path.c kthread.c $(LIBS) -lpthread


clean:
	rm -fr *.dSYM $(PROG)

