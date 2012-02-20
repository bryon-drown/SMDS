# $Id: Makefile,v 1.2 2006/05/18 22:06:45 mculbert Exp $

CFLAGS = -fPIC

random.so: random.o
	gcc -shared -o librandom.so random.o

clean :
	rm -f *~ librandom.so random.o smds/*~ smds/Cores/*~
