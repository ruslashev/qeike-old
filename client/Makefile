OBJS = $(patsubst ./%.cc,.objs/%.o, $(shell find . -type f -name '*.cc'))
warnings = -Wall -Wextra -Wshadow -Wno-unused-parameter -Wno-unused-variable \
		   -Wduplicated-cond -Wdouble-promotion -Wnull-dereference -Wformat=2 \
		   -Wdisabled-optimization -ftabstop=2 -Wsuggest-override -Wlogical-op
flags = -ggdb3 -Og -std=c++0x -fno-rtti
libraries = -lSDL2 -lGLEW -lGL
CXX = g++
CXXFLAGS = $(warnings) $(flags)
LDFLAGS = $(libraries)

all: objdir qeike
	./qeike

.objs/%.o: %.cc
	@echo "Compiling $<"
	@$(CXX) -c -o $@ $< $(CXXFLAGS)

qeike: $(OBJS)
	@echo "Linking to $@"
	@$(CXX) -o $@ $^ $(LDFLAGS)

valgrind: objdir qeike
	valgrind --leak-check=full ./qeike

callgrind: objdir qeike
	@valgrind --tool=callgrind ./qeike
	@kcachegrind callgrind.out.$!

objdir:
	@mkdir -p .objs

clean:
	@rm -f qeike
	@rm -f $(OBJS)

