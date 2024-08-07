TARGET = s21_matrix_oop.a
CC = gcc
CFLAGS = -Wall -Werror -Wextra -std=c++17 -O3 -lstdc++
TEST_FLAGS = -lgtest -pthread
TEST_TARGET = testing_exe
MODULES = $(wildcard *.cc)
OBJECTS = $(patsubst %.cc, %.o, $(MODULES))

all: clean $(TARGET)

build: $(TARGET)

rebuild: clean $(TARGET)

%.o: %.cc
	@mkdir -p obj
	@$(CC) -c $(CFLAGS) $< -o ./obj/$@

$(TARGET): $(OBJECTS)
	@ar rcs $(TARGET) $(wildcard obj/*.o)

test: $(TARGET)
	@g++ $(CFLAGS) ./tests/*.cc $(TEST_FLAGS) $(TARGET) -o ./tests/$(TEST_TARGET)
	@./tests/$(TEST_TARGET)

style_check:
	@cp ../materials/linters/.clang-format .
	@clang-format -n *.cc *.h **/*.cc **/*.h
	@rm ./.clang-format

format_style:
	@cp ../materials/linters/.clang-format .
	@clang-format -i *.cc *.h **/*.cc **/*.h
	@rm ./.clang-format

valgrind: test
	@CK_FORK=no valgrind --vgdb=no --leak-check=full \
	--show-leak-kinds=all --track-origins=yes --log-file="valgrind.log" -v --verbose -q --quiet -s tests/$(TEST_TARGET)

clean:
	@echo "Deleting unnecessary files..."
	@rm -rf obj *.a *.o tests/$(TEST_TARGET) *.dSYM **/*.dSYM *.log **/*.log

