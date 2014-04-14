CC = g++ -Wall -g

mmhc: Matrix.o main.o
	$(CC) -o $@ $+

Matrix.o: Matrix.cpp classes.h
	$(CC) -c -o $@ $<

main.o: main.cpp classes.h
	$(CC) -c -o $@ $<

clean:
	rm *.o