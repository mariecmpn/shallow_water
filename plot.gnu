# Graphique pour le modele de Shallow-Water

set title "Solution à T_{fin} pour les équations de Shallow-Water"
#plot 'solution_u.dat' with lines
plot 'solution_h.dat' with lines, 'topo.dat' with lines
#plot 'topo.dat' with lines