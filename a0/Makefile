#!/usr/bin/make

VECLIBDIR = $(PWD)/../vecmath
INCFLAGS  = -I /usr/include/GL
INCFLAGS += -I $(VECLIBDIR)/include

LINKFLAGS  = -lglut -lGL -lGLU
LINKFLAGS += -L $(VECLIBDIR) -lvecmath

CFLAGS    = -O2
CC        = g++
SRCS      = main.cpp
OBJS      = $(SRCS:.cpp=.o)
PROG      = a0

all: $(SRCS) $(PROG)

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(INCFLAGS) $(LINKFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -c -o $@ $(INCFLAGS)

depend:
	makedepend $(INCFLAGS) -Y $(SRCS)

clean:
	rm $(OBJS) $(PROG)

