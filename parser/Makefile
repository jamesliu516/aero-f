LIB = libparser.a

CCFLAGS = $(CXXFLAGS)

OBJS = Lexer.o \
       Parser.tab.o \
       Dictionary.o \
       AssignerCore.o

$(LIB): $(OBJS)
	$(SAR) $(LIB) $(OBJS)
	echo $(LIB)

Lexer.o: Parser.tab.c Lexer.C
	$(CXX) $(CXXFLAGS) -c Lexer.C

Parser.tab.o: Parser.tab.c
	$(CXX) $(CXXFLAGS) -c Parser.tab.c

Parser.tab.c: Parser.y
	bison -d -v -p yyCmdf Parser.y

Lexer.C: Lexer.l
	flex -oLexer.C -PyyCmdf Lexer.l

clean:
	rm -f $(LIB) $(OBJS) Lexer.C Parser.output Parser.tab.c Parser.tab.h
