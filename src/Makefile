#
# Makefile to build sequential and distributed versions of "lowx".
#

BINDIR  = $(HOME)/.local/bin

LEX     = flex
LEXLIB  = -lfl
YACC    = bison

CC      = cc
CFLAGS  = -I. -Wall -O3
LDFLAGS =
LIBS    = $(LEXLIB)

OBJS    = lowx.o \
          rcpr.o \
          cosetenum.o \
          tests.o \
          lowindex.o \
          conjugate.o \
          checkpoint.o \
          print.o \
          gens_init.o \
          lex.yy.o \
          parser_aux.o \
          parser.tab.o \
          print_error.o \

ifdef DEBUG
CFLAGS += -DDEBUG
endif

ifdef DIST
CFLAGS += -DDISTRIBUTED
OBJS   += distributed.o
endif

default:
	@echo Targets available:
	@echo
	@echo '  lowx          - Build sequential lowx.'
	@echo '  DIST=1 lowx   - Build the distributed lowx.'
	@echo '  DEBUG=1 lowx  - Build with debug mode on.'
	@echo "  install       - Install lowx into $(BINDIR)"
	@echo "  uninstall     - Delete $(BINDIR)/lowx"
	@echo '  tags          - Run ctags *.c *.h parser.y lexer.l'
	@echo '  clean         - Delete building artifacts.'
	@echo '  clobber       - clean + delete lowx executable.'
	@echo

lowx: $(OBJS)
	$(CC) $(LDFLAGS) -o lowx $(OBJS) $(LIBS)

%.o: %.c lowx.h
	$(CC) $(CFLAGS) -c $<

lowx.o: lowx.c lowx.h
	$(CC) $(CONFIG) $(CFLAGS) -c lowx.c

parser.tab.c parser.tab.h: parser.y lowx.h
	$(YACC) -d parser.y

parser.tab.o: parser.tab.c
	$(CC) $(CFLAGS) -c parser.tab.c

lex.yy.c: lexer.l
	$(LEX) lexer.l

lex.yy.o: lex.yy.c parser.tab.h
	$(CC) $(CFLAGS) -c lex.yy.c

install: $(BINDIR)/lowx

$(BINDIR)/lowx: lowx
	mkdir -p $(BINDIR)
	install lowx $(BINDIR)

uninstall:
	rm -f $(BINDIR)/lowx

clean:
	rm -f *.o *.tab.* lex.yy.c *.bak *~ core

clobber: clean
	rm -f lowx tags

