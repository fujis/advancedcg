# Compiler
COMPILER = g++
CXXFLAGS = -O3 -std=c++11

# library-related parameters
LDFLAGS = -L/opt/local/lib -lglfw -lGLEW -framework OpenGL
LIBS    = 

# Include folder specification
INCLUDE = -I../include -I./imgui

# Generated binary file name & directory
TARGETS = acg
TARGETDIR = ../bin

# Location of source files (it scan subfolders to make a list of *.cpp files)
SRCROOT   = .
SRCDIRS  := $(shell find $(SRCROOT) -type d)
SOURCES   = $(foreach dir, $(SRCDIRS), $(wildcard $(dir)/*.cpp))

# Location to put intermediate file (*.o) & its file name (it's automatically determined from cpp file)
OBJROOT   = .
OBJECTS   = $(addprefix $(OBJROOT)/, $(SOURCES:.cpp=.o)) 
OBJDIRS   = $(addprefix $(OBJROOT)/, $(SRCDIRS)) 

# Command to build a executable file
$(TARGETS): $(OBJECTS) $(LIBS)
	$(COMPILER) -o $(TARGETDIR)/$@ $^ $(LDFLAGS)

# Scannning the directory of intermediate files and compile *.cpp to intermediate file *.o, 
$(OBJROOT)/%.o: $(SRCROOT)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(COMPILER) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

run: $(TARGETS)
	cd $(TARGETDIR); ./$(TARGETS); cd -

clean:
	rm -f $(OBJECTS) $(TARGETDIR)/$(TARGETS)