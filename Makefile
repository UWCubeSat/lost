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
DOCUMENTATION := $(wildcard documentation/*.man)
DOCUMENTATION_FILES := $(patsubst documentation/man-%.h, documentation/%.man,$(DOCUMENTATION))
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

documentation/man-%.h: documentation/%.man
	xxd -i $< > $@

src/main.o: $(DOCUMENTATION_FILES)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -MMD -c $< -o $@

-include $(DEPS)

test: $(TEST_BIN)
	$(TEST_BIN)

$(TEST_BIN): $(TEST_OBJS)
	$(CXX) $(LDFLAGS) -o $(TEST_BIN) $(TEST_OBJS) $(LIBS)

clean:
	rm -f $(OBJS) $(DEPS) $(TEST_OBJS)
	rm -i $(BSC)

.PHONY: all clean test
