CC = gcc

INCDIR = include
OBJDIR = obj
# CFLAGS contains options passed to the compiler
CFLAGS = -c -I$(INCDIR)
OBJECTS = main.o $(OBJDIR)/evolution.o $(OBJDIR)/read_write_pgm_image.o $(OBJDIR)/utilities.o

gol.x: $(OBJECTS)
	$(CC) $(OBJECTS) -o $@
main.o: main.c
	$(CC) $(CFLAGS) main.c
$(OBJDIR)/%.o: src/%.c
	$(CC) $(CFLAGS) $^ -o $@
clean:
	rm -rf *.o *.x
	rm -rf $(OBJDIR)/*.o
