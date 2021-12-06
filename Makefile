CC = g++
FLAG = -fpermissive
SOURCE = src/
MAIN = main.out

build: 
	@echo "Compile" 
	${CC} ${FLAG} ${SOURCE}*.cpp -o src/${MAIN}

run:
	@echo "Run"
	src/./${MAIN}

clean:
	@echo "Clean up this shit man"
	rm src/*.out