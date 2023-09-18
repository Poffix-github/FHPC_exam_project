CC = mpicc

INCDIR = include
OBJDIR = obj
# CFLAGS contains options passed to the compiler
CFLAGS = -c -openmp -I$(INCDIR)
OBJECTS = main.o $(OBJDIR)/evolution_parallel.o $(OBJDIR)/read_write_pgm_image.o $(OBJDIR)/utilities.o

gol.x: $(OBJECTS)
	$(CC) $(OBJECTS) -o $@
main.o: main.c
	$(CC) $(CFLAGS) main.c -o $@
$(OBJDIR)/%.o: src/%.c
	$(CC) $(CFLAGS) $^ -o $@
clean:
	rm -rf *.o *.x
	rm -rf $(OBJDIR)/*.o
