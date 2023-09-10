CC=g++
CFLAGS=-c -v -mtune=generic
SOURCES=main.cpp enter.cpp calculate_answer.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS)