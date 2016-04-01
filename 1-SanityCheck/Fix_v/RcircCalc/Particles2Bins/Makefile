CC = gcc

OPTIMIZE = -Wall -O3

OPTIONS = $(OPTIMIZE)
NAME = particles2bins
EXEC = $(NAME).x

OBJS = $(NAME).o


INCL = readgadget.h 

CFLAGS = $(OPTIONS)

LIBS = -lm

$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)

$(OBJS):$(INCL)

clean:
	rm -f $(OBJS) $(EXEC)
