CC = gcc

all: t_trnpp

t_trnpp: trnpp.o main_trnpp.o
	$(CC) -o t_trnpp trnpp.o main_trnpp.o

trnpp.o: trnpp.c
	$(CC) -c trnpp.c
main_trnpp.o: main_trnpp.c
	$(CC) -c main_trnpp.c

clean:
	-rm trnpp.o
	-rm main_trnpp.o
