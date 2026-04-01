# Compiler
CXX = g++
CXXFLAGS = -std=c++17 -Iinclude

# Source files
SRC = src/matrix.cpp src/gaussian_solver.cpp src/lu_solver.cpp src/system_linear_eq.cpp  src/iterative_solver.cpp

# Object files
OBJ = $(SRC:.cpp=.o)

# Static and shared libraries
STATIC_LIB = libmatrix.a
SHARED_LIB = libmatrix.so

# Default target
all: $(STATIC_LIB) $(SHARED_LIB) matrix_static matrix_dynamic

# Compile object files
src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -fPIC -c $< -o $@

# Create static library
$(STATIC_LIB): $(OBJ)
	ar rcs $(STATIC_LIB) $(OBJ)

# Create shared library
$(SHARED_LIB): $(OBJ)
	$(CXX) -shared -o $(SHARED_LIB) $(OBJ)

# Build executable with static library
matrix_static: main.cpp $(STATIC_LIB)
	$(CXX) $(CXXFLAGS) main.cpp libmatrix.a -o matrix_static

# Build executable with dynamic library
matrix_dynamic: main.cpp $(SHARED_LIB)
	$(CXX) $(CXXFLAGS) main.cpp -L. -lmatrix -o matrix_dynamic

# Run program
run_static:
	./matrix_static

run_dynamic:
	LD_LIBRARY_PATH=. ./matrix_dynamic

# Clean build files
clean:
	rm -f src/*.o $(STATIC_LIB) $(SHARED_LIB) matrix_static matrix_dynamic