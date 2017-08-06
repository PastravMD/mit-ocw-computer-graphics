#!/usr/bin/make

TARGET = vecmath.so
LIBDIR = $(PWD)
SOURCES = $(shell echo $(LIBDIR)/src/*.cpp)
INCLUDE = $(LIBDIR)/include/
OBJECTS = $(SOURCES:.cpp=.o)

CC = gcc
CC_FLAGS = -c -Wall -fPIC -O3 -I$(INCLUDE)
LD_FLAGS = -shared

.PHONY : all 
all: $(TARGET)

$(TARGET) : $(OBJECTS)
	$(CC) $(LD_FLAGS) $(OBJECTS) -o $@ 
	ln -s $@ lib$@

$(SOURCES:.cpp=.o):%o:%cpp
	$(CC) $(CC_FLAGS) $< -o $@

.PHONY : clean
clean:
	rm -f $(OBJECTS) $(TARGET)
