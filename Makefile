# Copyright (c) 2020 Mark Polyakov (If you edit the file, add your name here!)
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Simple makefile: Compile all .c files into .o files, generating "dependency" .d files too (see
# https://stackoverflow.com/q/2394609)

SRCS := $(wildcard src/*.cpp)
TESTS := $(wildcard test/*.cpp)
MANS := $(wildcard documentation/*.man)
MAN_TXTS := $(patsubst documentation/%.man, documentation/%.txt, $(MANS))
MAN_HS := $(patsubst documentation/%.man, documentation/man-%.h, $(MANS))
OBJS := $(patsubst %.cpp,%.o,$(SRCS))
TEST_OBJS := $(patsubst %.cpp,%.o,$(TESTS) $(filter-out %/main.o, $(OBJS)))
DEPS := $(patsubst %.cpp,%.d,$(SRCS) $(TESTS)) # includes tests
BIN  := lost
TEST_BIN := ./lost-test

BSC  := bright-star-catalog.tsv

LIBS     := -lcairo
CXXFLAGS := $(CXXFLAGS) -Ivendor -Isrc -Idocumentation -Wall -pedantic --std=c++11

all: $(BIN) $(BSC)

$(BSC): download-bsc.sh
	./download-bsc.sh

$(BIN): $(OBJS)
	$(CXX) $(LDFLAGS) -o $(BIN) $(OBJS) $(LIBS)

documentation/%.txt: documentation/%.man
	groff -mandoc -Tascii $< > $@

documentation/man-%.h: documentation/%.txt
	xxd -i $< > $@

src/main.o: $(MAN_HS)

docs:
	doxygen

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -MMD -c $< -o $@

-include $(DEPS)

test: $(TEST_BIN)
	$(TEST_BIN)
	make script-tests

$(TEST_BIN): $(TEST_OBJS)
	$(CXX) $(LDFLAGS) -o $(TEST_BIN) $(TEST_OBJS) $(LIBS)

script-tests:
	bash ./test/scripts/pyramid-incorrect.sh
	bash ./test/scripts/readme-examples-test.sh
	rm -f img_7660.png my-database.dat attitude.txt annotated-7660.png annotated-input.png raw-input.png input.png img_7660.png.1

clean:
	rm -f $(OBJS) $(DEPS) $(TEST_OBJS) $(MAN_HS)
	rm -f img_7660.png my-database.dat attitude.txt annotated-7660.png annotated-input.png raw-input.png input.png img_7660.png.1

clean_all: clean
	rm -f $(BSC)

.PHONY: all clean test
