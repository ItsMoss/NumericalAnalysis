CC = g++
CFLAGS = -Wall -Werror -std=gnu++98 -pedantic -ggdb3

numerics: funky.o helpers.o numerics.o
	$(CC) $(CFLAGS) funky.o helpers.o numerics.o -o numerics

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

clean:
	rm -rf *.o numerics
