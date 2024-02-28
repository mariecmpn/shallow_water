# Graphique pour le modele de Shallow-Water

set title "Solution à T_{fin} pour les équations de Shallow-Water"
#plot 'solution_u.dat' with lines lw 2 linecolor rgb "light-magenta" title "vitesse"
plot 'solution_h.dat' with lines lw 2 linecolor rgb "light-turquoise" title "hauteur d'eau", 'topo.dat' with lines lw 2 linecolor rgb "forest-green" title "topographie"
#plot 'solution_h.dat' with lines lw 2 linecolor rgb "light-turquoise" title "hauteur d'eau"