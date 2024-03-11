# Graphique pour le modele de Shallow-Water

set title "Solution à T_{fin} pour les équations de Shallow-Water"
plot 'solutions.dat' u 1:3 with lines lw 2 linecolor rgb "light-magenta" title "vitesse"
#plot 'solutions.dat' u 1:2 with lines lw 2 linecolor rgb "light-turquoise" title "hauteur d'eau", 'topo.dat' with lines lw 2 linecolor rgb "forest-green" title "topographie"
#plot 'solutions.dat' u 1:4 with lines lw 2 linecolor rgb "plum" title "débit"