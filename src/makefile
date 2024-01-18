# Compilateur utilise
FC = gfortran-11


# Options de compilation
FFLAGS = -O0

# dossier pour les fichiers binaires
ODIR = ../bin

# chemin pour les fichiers executables
EXEDIR = ..

# Liste des fichiers binaires
OBJS = $(ODIR)/main.o $(ODIR)/schemas.o $(ODIR)/initialisation.o $(ODIR)/numerics.o

# Pour le modele scalaire
shallowwater: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $(EXEDIR)/trafic_scalaire


# regle clean pour supprimer les fichiers binaires et executable
clean:
	rm $(ODIR)/*.o $(EXEDIR)/shallowwater

# Regles construction des fichiers binaires
$(ODIR)/main.o: main.f90 schemas.f90 init_posttraitement.f90 numerics.f90 
	$(FC) -c $(FFLAGS) main.f90 -o $(ODIR)/main.o
$(ODIR)/numerics.o: numerics.f90 
	$(FC) -c $(FFLAGS) numerics.f90 -o $(ODIR)/numerics.o
$(ODIR)/schemas.o: schemas.f90 init_posttraitement.f90 numerics.f90 
	$(FC) -c $(FFLAGS) schemas.f90 -o $(ODIR)/schemas.o
$(ODIR)/init_posttraitement.o: init_posttraitement.f90 numerics.f90 
	$(FC) -c $(FFLAGS) init_posttraitement.f90 -o $(ODIR)/init_posttraitement.o