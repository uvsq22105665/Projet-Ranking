#graphe et méthode à tester : puissances ou aitken
GRAPHE = Stanford.txt
METHODE = puissances


DEBUG = no


	CFLAGS =


#Nom de l'exécutable
EXEC = pagerank

#Liste des fichiers sources et objets
SRC = $(wildcard *.c)
OBJ = $(SRC:.o=.c)


#Edition de liens / génération de l'exécutable
all: $(EXEC)

	./$(EXEC) $(METHODE) graphes/$(GRAPHE) > graphes/$(METHODE)-$(GRAPHE)


$(EXEC) : $(OBJ)
	@gcc -o $@ $^ -lm

#Génération des graphes
plot:
	gnuplot normes.p

#Compilation des .c
main.o : matrice.h puissances.h aitkenprim.h

%.o : %.c def.h
	@gcc $(CFLAGS) -o $@ -c $< 
	

#Nettoyage du projet
.PHONY : $(EXEC) plot clean mrproper

clean:
	rm -rf *.o

mrproper: clean
	rm -rf pagerank
