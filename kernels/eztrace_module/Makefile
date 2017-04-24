# Copyright (C) CNRS, INRIA, Université Bordeaux 1, Télécom SudParis
# See COPYING in top-level directory.

prefix=/home/pichon/eztrace-1.1-5_install
exec_prefix=${prefix}
libdir=${exec_prefix}/lib
includedir=${prefix}/include

BIN = libeztrace-convert-kernels.so libeztrace-kernels.so libeztrace-autostart-kernels.so

ldflags=$(LDFLAGS)  -L${libdir}
cflags=$(CFLAGS) -I/home/pichon/eztrace-1.1-5_install/include -DGTG_OUT_OF_ORDER -I/home/pichon/eztrace-1.1-5_install/include/ -DTID_RECORDING_ENABLED -I${includedir}

CONVERT_LDFLAGS=-L/home/pichon/eztrace-1.1-5_install/lib -L/home/pichon/eztrace-1.1-5_install/lib/
RECORD_LDFLAGS=-L/home/pichon/eztrace-1.1-5_install/lib/

all: $(BIN)

libeztrace-convert-kernels.so: eztrace_convert_kernels.c kernels_ev_codes.h
	$(CC) --shared -o libeztrace-convert-kernels.so eztrace_convert_kernels.c -fPIC $(cflags) $(ldflags) $(CONVERT_LDFLAGS) -ldl

libeztrace-kernels.so: /home/pichon/pastix_eztrace/build/kernels/libpastix_kernels.a kernels_ev_codes.h
	$(CC) --shared -o libeztrace-kernels.so -fPIC $(cflags) $(ldflags) -Wl,-Bsymbolic $(RECORD_LDFLAGS) -ldl /home/pichon/pastix_eztrace/build/kernels/libpastix_kernels.a

libeztrace-autostart-kernels.so: /home/pichon/pastix_eztrace/build/kernels/libpastix_kernels.a kernels_ev_codes.h
	$(CC) --shared -o libeztrace-autostart-kernels.so -fPIC $(cflags) $(ldflags) -Wl,-Bsymbolic $(RECORD_LDFLAGS) -ldl -DEZTRACE_AUTOSTART /home/pichon/pastix_eztrace/build/kernels/libpastix_kernels.a

clean:
	rm -f $(BIN)
