CC ?= gcc
CFLAGS ?= -O2 -Wall -Wextra -std=c11
PKG_LIBS := $(shell pkg-config --libs msolve)
PKG_CFLAGS := $(shell pkg-config --cflags msolve)

msolve-inproc: msolve_inproc.c
	$(CC) $(CFLAGS) $(PKG_CFLAGS) -o $@ $< $(PKG_LIBS)

clean:
	rm -f msolve-inproc

.PHONY: clean
