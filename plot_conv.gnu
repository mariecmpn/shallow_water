# Graphique de convergence pour Shallow-Water

set logscale y
set logscale x
set xlabel "Δx"
set ylabel "erreur"
set title "Graphique de convergence pour ||h_i^n+Z_i-H||"
plot "erreursL2.dat" using 1:2 with lines lw 2 title "erreur L^2", "erreursinf.dat" using 1:2 with lines lw 2 title "erreur L^∞", "erreursL2.dat" using 1:1 with lines lw 2 title "O(Δx)"
#set title "Graphique de convergence pour ||u_i^n||"
#plot "erreursL2.dat" using 1:3 with lines lw 2 title "erreur L^2", "erreursinf.dat" using 1:3 with lines lw 2 title "erreur L^∞", "erreursL2.dat" using 1:1 with lines lw 2 title "O(Δx)"