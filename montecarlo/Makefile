.PHONY: all clean
CC = gcc
LD = gcc
LD_FLAGS = -shared
C_FLAGS = -g -Ofast -fPIC -std=gnu99


all: pressure.so fermi_gas.e lj_gas.e NM/NM.e NM/layers/NM_layers.e NSM/NSM_layers.e fermi_gas_red.e params_Pauli/fermi_gas_params.e

%.e: %.o
	$(CC) $(C_FLAGS) $^ -o $@ -lm

%.so: %.o
	$(LD) $(LD_FLAGS) $^ -o $@

%.o: %.c
	$(CC) $(C_FLAGS) $^ -c -o $@

clean:
	rm -rfv lj.so lj.o morse.so morse.o NM/NM.o
