# Source and build directories
SRCDIR = src
OBJDIR = obj

# Source and object files
SRCS = $(wildcard $(SRCDIR)/*.c)
OBJS = $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(SRCS))

# Target executable name
TARGET = sw.out

# Compiler and flags
CC = gcc
CFLAGS = -O3 -march=native -Wall -Wextra -fopenmp -std=c11 -Isrc -lm
LDFLAGS = -fopenmp

# Libraries
LDLIBS = -lm


# Default target
all: $(TARGET)

# Clean up build files
clean:
	rm -rf $(OBJDIR) $(TARGET)

# Rebuild everything from scratch
rebuild: clean all

# Phony targets
.PHONY: all clean rebuild


# Rule to build the target executable
$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS) $(LDLIBS)

# Rule to build object files
$(OBJDIR)/%.o: $(SRCDIR)/%.c | $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Ensure the object directory exists
$(OBJDIR):
	mkdir -p $(OBJDIR)

