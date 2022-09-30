CFLAGS=		-g3 -O2 -Wall -Wno-unused -Wno-unused-result -Wno-sign-compare
CC=             gcc
CPPFLAGS=
INCLUDES=	
OBJS=		yak-priv.o count.o compare.o htab.o kthread.o 
PROG=		yam
LIBS=		-lm -lz -lpthread

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.SUFFIXES:.c .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

yam:$(OBJS) main.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o ext/*.o a.out $(PROG) *~ *.a *.dSYM session*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE



htab.o: kthread.h yak-priv.h khashl.h ksort.h
main.o: ketopt.h yak-priv.h count.h
count.o: kthread.h yak-priv.h count.h kseq.h
compare.o: khashl.h htab.h yak-priv.h
kthread.o: kthread.h
