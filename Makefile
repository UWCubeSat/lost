# Copyright (c) 2020 Mark Polyakov, Karen Haining, Muki Kiboigo (If you edit the file, add your name here!)
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
DOXYGEN_DIR := ./documentation/doxygen
OBJS := $(patsubst %.cpp,%.o,$(SRCS))
TEST_OBJS := $(patsubst %.cpp,%.o,$(TESTS) $(filter-out %/main.o, $(OBJS)))
DEPS := $(patsubst %.cpp,%.d,$(SRCS) $(TESTS)) # includes tests
BIN  := lost
TEST_BIN := ./lost-test

BSC  := bright-star-catalog.tsv

LIBS     := -lcairo
CXXFLAGS := $(CXXFLAGS) -Ivendor -Isrc -Idocumentation -Wall -Wextra -Wno-missing-field-initializers -pedantic --std=c++11
RELEASE_CXXFLAGS := $(CXXFLAGS) -O3
# debug flags:
CXXFLAGS := $(CXXFLAGS) -ggdb -fno-omit-frame-pointer
ifndef LOST_DISABLE_ASAN
	CXXFLAGS := $(CXXFLAGS) -fsanitize=address
endif

RELEASE_LDFLAGS := $(LDFLAGS)

# debug link flags:
ifndef LOST_DISABLE_ASAN
	LDFLAGS := $(LDFLAGS) -fsanitize=address
endif

# Use Double Mode by default.
# If compiled with LOST_FLOAT_MODE=1, we will use floats.
ifdef LOST_FLOAT_MODE
	CXXFLAGS := $(CXXFLAGS) -Wdouble-promotion -Werror=double-promotion -D LOST_FLOAT_MODE
endif

all: $(BIN) $(BSC)

release: CXXFLAGS := $(RELEASE_CXXFLAGS)
release: LDFLAGS := $(RELEASE_LDFLAGS)
release: all

$(BIN): $(OBJS)
	$(CXX) $(LDFLAGS) -o $(BIN) $(OBJS) $(LIBS)

documentation/%.txt: documentation/%.man
	groff -mandoc -Tascii $< > $@
	printf '\0' >> $@

documentation/man-%.h: documentation/%.txt
	xxd -i $< > $@

src/main.o: $(MAN_HS)

docs:
	doxygen

lint:
	cpplint --recursive src test

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -MMD -c $< -o $@

-include $(DEPS)

test: $(BIN) $(BSC) $(TEST_BIN)
	$(TEST_BIN)
	# bash ./test/scripts/pyramid-incorrect.sh
	bash ./test/scripts/readme-examples-test.sh
	bash ./test/scripts/random-crap.sh

$(TEST_BIN): $(TEST_OBJS)
	$(CXX) $(LDFLAGS) -o $(TEST_BIN) $(TEST_OBJS) $(LIBS)

clean:
	rm -f $(OBJS) $(DEPS) $(TEST_OBJS) $(MAN_HS)
	rm -rf $(DOXYGEN_DIR)

clean_all: clean
	rm -f $(BSC)

.PHONY: all clean test docs lint
