CC = gcc
CFLAGS = -O3
LDLIBS = -lm -fopenmp -Wno-alloc-size-larger-than -g
RM = rm -f
OBJS = project.o
EXECUTABLE = project
all:main.o

execute:
		./$(EXECUTABLE)
main.o:
		$(CC) -o project project.c $(LDLIBS) $(CFLAGS)

clean:
		$(RM) $(OBJS) $(EXECUTABLE)