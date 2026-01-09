TARGETS = $(patsubst %.c,lib%.so,$(SOURCES))
SOURCES = $(wildcard *.c)

PYTHONCFLAGS  = $(shell python3-config --cflags)
PYTHONLDFLAGS = $(shell python3-config --ldflags) -shared

CFLAGS  = -I. -Wall -fPIC
LDFLAGS = -lm

CC = /usr/bin/gcc

# VERBOSE := @   # uncomment to hide commands
VERBOSE :=

all: $(TARGETS)

# Rule: libX.so from X.c
lib%.so: %.c
	$(VERBOSE)$(CC) $(CFLAGS) $(PYTHONCFLAGS) -shared -o $@ $< $(LDFLAGS) $(PYTHONLDFLAGS)

# Optional: keep your existing “binary from .c” rule
%: %.c
	$(VERBOSE)$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

clean:
	rm -f core a.out *~ *.o $(TARGETS)
	rm -rf .codelite .clang