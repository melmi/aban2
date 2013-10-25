SRCDIR = src
OBJDIR = bin
DEPDIR = bin
BINDIR = bin

CXX      = g++
CXXFLAGS = -Wall -O3 -std=c++11
LDFLAGS  = -ljsoncpp

TARGET = $(BINDIR)/aban2
SRCS   = $(wildcard $(SRCDIR)/*.cpp)
OBJS   = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCS))
DEPS   = $(patsubst $(SRCDIR)/%.cpp,$(DEPDIR)/%.depends,$(SRCS))

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(OBJS)
	@mkdir -p $(BINDIR)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LDFLAGS) -o $(TARGET)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp 
	@mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(DEPDIR)/%.depends: $(SRCDIR)/%.cpp 
	@mkdir -p $(DEPDIR)
	$(CXX) -MM $(CXXFLAGS) $< -MT $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$<) -MF $@

clean:
	@rm -f $(OBJS) $(DEPS) $(TARGET)

run:
	$(TARGET)

-include $(DEPS)
