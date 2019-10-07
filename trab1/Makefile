    CC     = gcc -std=c11
    CFLAGS =
    LFLAGS = -lm

      PROG = pdeSolver
      OBJS = utils.o \
             edp.o \
	     gaussSeidel.o \
             $(PROG).o


.PHONY: limpa faxina clean distclean purge all

%.o: %.c %.h
	$(CC) -c $(CFLAGS) $(LFLAGS) $<

$(PROG):  $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)

clean:
	@rm -f *~ *.bak

purge:   clean
	@rm -f *.o core a.out
	@rm -f $(PROG)

