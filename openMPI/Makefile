IDIR = nqueue
CC = mpic++
CFLAGS = -I $(IDIR) -g -O3
SRC_DIR = .

# List all of the object files that need to be compiled
OBJ = tsp.o

# Specify the name of the executable
EXEC = tsp

# Build the executable by linking the object files together
$(EXEC): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^

# Compile each source file into an object file
%.o: $(SRC_DIR)/%.cpp $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

# Clean up object files and the executable
clean:
	rm -f $(OBJ) $(EXEC)
