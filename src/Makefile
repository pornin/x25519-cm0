CC = arm-linux-gcc
CFLAGS = -Wall -Wextra -Wshadow -Wundef -Os -mcpu=cortex-m0plus
LD = arm-linux-gcc
LDFLAGS =
LIBS =

OBJ = x25519-cm0.o
OBJ_TESTS = test_x25519.o

all: test_x25519

clean:
	-rm -f $(OBJ) $(OBJ_TESTS) test_x25519 test_x25519.gdb

test_x25519: $(OBJ) $(OBJ_TESTS)
	$(LD) $(LDFLAGS) -o test_x25519 $(OBJ) $(OBJ_TESTS) $(LIBS)

x25519-cm0.o: x25519-cm0.S
	$(CC) $(CFLAGS) -c -o x25519-cm0.o x25519-cm0.S

test_x25519.o: test_x25519.c
	$(CC) $(CFLAGS) -c -o test_x25519.o test_x25519.c
