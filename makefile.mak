CC=g++
CFLAGS=-I.
DEPS = frequency.h
OBJ = main.o frequency.o 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

frequency: $(OBJ)
	g++ -o $@ $^ $(CFLAGS) -lboost_iostreams