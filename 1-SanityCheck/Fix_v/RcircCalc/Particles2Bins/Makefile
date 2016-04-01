CC = gcc

OPTIMIZE = -Wall -O3

OPTIONS = $(OPTIMIZE)

EXEC = particles2bins

OBJS = ${EXEC}.o


INCL = readgadget.h 

CFLAGS = $(OPTIONS)

LIBS = -lm

$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)

$(OBJS):$(INCL)

clean:
	rm -f $(OBJS) $(EXEC)
