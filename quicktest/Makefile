CXX      = g++
CXXFLAGS = -std=c++17 -DDEBUG -g

SRCDIR = .
BINARY = test

%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $^

test: quicktest.o
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm *.o && rm mmu0 && rm mmu1