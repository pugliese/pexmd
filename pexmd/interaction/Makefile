.PHONY: all clean
CC = gcc
LD = gcc
LD_FLAGS = -shared
C_FLAGS = -g -Ofast -fPIC -std=gnu99


all: lj.so morse.so pauli.so coulomb.so qcnm.so

%.so: %.o
	$(LD) $(LD_FLAGS) $^ -o $@

%.o: %.c
	$(CC) $(C_FLAGS) $^ -c

clean:
	rm -rfv lj.so lj.o morse.so morse.o
