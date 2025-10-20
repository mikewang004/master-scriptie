
TARGETS= find_box_id.so 

PYTHONCFLAGS=$(shell python3-config --cflags)
PYTHONLDFLAGS=$(shell python3-config --cflags) -shared

CFLAGS=-I. -Wall -fPIC
LDFLAGS=-lm

CC=/usr/bin/gcc

# Use @ here to suppress all output :
#VERBOSE := @
 VERBOSE :=

all: $(TARGETS)

%.so : %.c
	$(VERBOSE)$(CC) $(CFLAGS) $(PYTHONCFLAGS)  $< -o $@ $(LDFLAGS) $(PYTHONLDFLAGS)

% : %.c
	$(VERBOSE)$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

$(OBJS) : Makefile

clean:
	rm -f core a.out *~ *.o $(TARGETS)
	rm -rf .codelite .clang