SRCDIR = src
OBJDIR = bin
DEPDIR = bin
BINDIR = bin
OUTDIR = out

CXX      = g++
CXXFLAGS = -O3 -std=c++11 -Wall
# CXXFLAGS = -g -O0 -std=c++11 -Wall
LDFLAGS  = -ljsoncpp -Xlinker -zmuldefs

TARGET = $(BINDIR)/aban2
SRCS   = $(wildcard $(SRCDIR)/*.cpp)
OBJS   = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCS))
DEPS   = $(patsubst $(SRCDIR)/%.cpp,$(DEPDIR)/%.depends,$(SRCS))

.PHONY: all clean run cleanout

all: $(TARGET)
	@echo Done.

$(TARGET): $(OBJS)
	@echo Linking $@
	@mkdir -p $(BINDIR)
	@mkdir -p $(OUTDIR)
	@$(CXX) $(CXXFLAGS) $(OBJS) $(LDFLAGS) -o $(TARGET)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp 
	@echo Compiling $<
	@mkdir -p $(OBJDIR)
	@$(CXX) $(CXXFLAGS) -c $< -o $@

$(DEPDIR)/%.depends: $(SRCDIR)/%.cpp 
	@echo Generating dependencies $<
	@mkdir -p $(DEPDIR)
	@$(CXX) -MM $(CXXFLAGS) $< -MT $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$<) -MF $@

clean:
	@rm -f $(OBJS) $(DEPS) $(TARGET)

cleanout:
	@rm -f $(OUTDIR)/*

run: all
	$(TARGET)

-include $(DEPS)
