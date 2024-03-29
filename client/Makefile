warnings = -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable \
		   -Wduplicated-cond -Wdouble-promotion -Wnull-dereference -Wformat=2 \
		   -Wdisabled-optimization -Wsuggest-override -Wlogical-op
flags = -ggdb3 -Og -std=c++0x -fno-rtti
libraries = -lSDL2 -lGLEW -lGL
CXX = g++
BIN = qeike

OBJS = $(patsubst ../%.cc,.objs/%.o, $(shell find .. -type f -name '*.cc'))
DEPS = $(OBJS:.o=.d)
CXXFLAGS = $(warnings) $(flags)
LDFLAGS = $(libraries)

$(shell mkdir -p .objs >/dev/null)
$(shell mkdir -p .objs/engine >/dev/null)
$(shell mkdir -p .objs/client >/dev/null)

all: $(BIN)
	./$(BIN)

$(BIN): $(OBJS)
	@echo "Linking to $@"
	@$(CXX) -o $@ $^ $(LDFLAGS)

.objs/client/%.o: %.cc
	@echo "Compiling $<"
	@$(CXX) -MMD -MP -c -o $@ $< $(CXXFLAGS)

.objs/engine/%.o: ../engine/%.cc
	@echo "Compiling $<"
	@$(CXX) -MMD -MP -c -o $@ $< $(CXXFLAGS)

gdb: $(BIN)
	gdb $(BIN)

valgrind: $(BIN)
	valgrind --leak-check=full ./$(BIN)

callgrind: $(BIN)
	@valgrind --tool=callgrind ./$(BIN)
	@kcachegrind callgrind.out.$!

.PHONY : clean
clean:
	@rm -f $(BIN) $(OBJS) $(DEPS)

-include $(DEPS)

