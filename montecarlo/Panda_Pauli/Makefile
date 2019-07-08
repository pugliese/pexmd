.PHONY: default help object executable all clean
CC = gcc

LD = $(CC)

SOURCE_C = $(wildcard src/*.c)
OBJECTS_C = $(patsubst src/%.c, bin/%_c.o, $(SOURCE_C))


CFLAGS = -c -fPIC -std=gnu99 -Ofast -Wall
EXECUTABLE = NM_PP.e
LIBRARY = lib.so

default: help

help:
	@echo $(SOURCE_C)
	@echo $(OBJECTS_C)
	@echo "\
Options:\n\n\
  make objects:       compiler makes objects for every *.c\n\
  make executable:    compiler makes executable\n\
  make all:           build all previous\n\
  make clean:         delete output files\n\
  make help:          display this help"

objects: $(OBJECTS_C)

lib: $(LIBRARY)

executable: $(EXECUTABLE)

all: objects executable lib

bin/%_c.o: src/%.c
	$(CC) $(CFLAGS) -c $^ -o $@ -lm

%.e: $(OBJECTS_C)
	$(LD) $^ -o $@ -lm

%.so: $(OBJECTS_C)
	$(CC) -shared $^ -o $@

clean:
	rm -rfv $(OBJECTS_C) $(EXECUTABLE1)
