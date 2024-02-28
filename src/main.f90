program systeme_SW
    use numerics
    use initialisation_sauvegarde
    use schemasSW
    IMPLICIT NONE

    ! Definition des variables
    real(rp) :: x_deb, x_fin ! dimensions du domaine en x: [x_deb,x_fin]
    integer :: Ns, Ns1, nb, dn ! nombre de cellules, nombre max de cellules si courbe de convergence, nb de nombres de cellules a tester si courbe conv, delta_n entre chaque nombre de cellules a tester different
    real(rp) :: CFL, T_fin, date, dt, dx ! condition CFL, Temps final, date pour laquelle on calcule l'iteration, delta_t, delta_x
    real(rp) :: v, Delta ! vitesse pour calculer dt avec la CFL, Delta = dt/dx pour ne pas avoir a le calculer tout le temps
    integer :: i,j ! entiers pour les boucles do
    real(rp) :: uL, uR, hL, hR ! conditions initiales
    real(rp), dimension(:,:), allocatable :: W_O ! tableau W^n des solutions
    real(rp), dimension(:,:), allocatable :: W_Op, W_Om ! tableaux pour la reconstruction hydrostatique
    real(rp), dimension(:,:), allocatable :: W_N ! tableau W^{n+1} des solutions
    real(rp), dimension(:,:), allocatable :: Flux ! tableau des flux numeriques
    real(rp), dimension(:), allocatable :: Zi ! tableau pour la topographie
    real(rp), dimension(:), allocatable :: Err_u, Err_h ! tableaux pour les erreurs (si pas de graphe de convergence)
    real(rp), dimension(:,:), allocatable :: Err ! tableaux pour les erreurs (si graphe de convergence)
    integer, dimension(:), allocatable :: N ! tableaux pour les Ns si graphe de convergence
    character(len = 1) :: conv ! convergence ou non
    real(rp), dimension(2,2) :: cond
    character(len = 2) :: schema ! schema utilise
    integer :: Nb_iter = 0 ! nombre d'iterations en temps
    integer :: topo ! quelle topographie on utilise

    write(6,*) '------------------------------------------'
    write(6,*) '----------- Resolution syteme ------------'
    write(6,*) '------------- SHALLOW-WATER --------------'
    write(6,*) '------------------------------------------'

    write(6,*)

    ! lecture des donnees du fichier donnees.dat
    call lecture_donnees_syst('init.dat',x_deb,x_fin,Ns,CFL,T_fin,cond,&
                                schema,uL,uR,hL,hR,topo,conv,Ns1,nb)

    write(6,*)

    if (conv == 'Y') then ! si on veut faire des graphiques de convergence
        write(6,*) 'CALCUL GRAPHIQUES DE CONVERGENCE'
        write(6,*) 'On calcule les erreurs pour n nombres de cellules différents entre Ns_deb et Ns_fin, avec: '
        write(6,*) 'n = ', nb
        write(6,*) 'Ns_deb = ', Ns
        write(6,*) 'Ns_fin = ', Ns1
        write(6,*)

        allocate(N(nb),Err(3,nb))
        dn = (Ns1 - Ns)/(nb-1)
        do i = 1,nb
            N(i) = Ns + (i-1)*dn
        end do

        do j = 1,nb
            write(6,*) 'Calcul pour Ns = ', N(j)
            ! calcul pas de maillage en espace dx
            dx = (x_fin - x_deb)/N(j)

            ! allocation memoire des tableaux
            allocate(W_O(2,N(j)), W_N(2,N(j)), Flux(2,N(j)-1), Err_u(N(j)), Err_h(N(j)), Zi(N(j)))
            allocate(W_Om(2,N(j)),W_Op(2,N(j)))

            ! calcul topographie
            call topographie(Zi, N(j), dx, x_deb, topo)

            ! initialisation pour t = 0
            call initialisation_syst(W_O, N(j), x_deb, x_fin, uL, uR, hL, hR, Zi) ! en variables primitives

            ! boucle en temps
            date = 0._rp
            do while (date < T_fin)
                ! CFL 
                v = abs(lambda_1(W_O(:,1))) ! on fait le max des valeurs propres
                v = max(v,abs(lambda_2(W_O(:,1))))
                do i = 2,N(j)
                    v = max(v, abs(lambda_1(W_O(:,i))), abs(lambda_2(W_O(:,i))))
                end do
                dt = dx*CFL/v
                dt = min(dt, T_fin - date)
                date = date + dt

                ! on etait en variables primitives donc on passe en variables conservatives
                call prim_to_conserv(W_O, N(j))
                
                ! calcul des flux
                if (schema == 'LF') then
                    call flux_LF_syst(N(j), Flux, W_O, dt, dx)
                else if (schema == 'RS') then
                    call flux_RS_syst(N(j), Flux, W_O)
                else if (schema == 'HL') then
                    call flux_HLL_syst(N(j), Flux, W_O, dx, dt)
                else if (schema == 'HY') then
                    call flux_recons_hydro(N(j), Flux, W_O, Zi, dx, dt, W_Om, W_Op)
                end if

                ! update calcul de u_i^{n+1}
                Delta = dt/dx
                if (schema == 'HY') then
                    do i = 2,(N(j)-1)
                        W_N(1,i) = W_O(1,i) - Delta*(Flux(1,i) - Flux(1,i-1))
                        W_N(2,i) = W_O(2,i) - Delta*(Flux(2,i) - Flux(2,i-1)) + dt*terme_src_hy(dx,W_Op(1,i-1),W_Om(1,i))
                    end do
                else
                    do i = 2,(N(j)-1)
                        W_N(1,i) = W_O(1,i) - Delta*(Flux(1,i) - Flux(1,i-1))
                        W_N(2,i) = W_O(2,i) - Delta*(Flux(2,i) - Flux(2,i-1)) + dt*terme_src(W_O(1,i),dx,Zi(i-1), Zi(i+1))
                    end do
                end if

                ! Conditions aux limites
                !if (condition == 'D') then 
                !    ! Dirichlet
                !    W_N(1,1) = W_O(1,1)
                !    W_N(1,N(j)) = W_O(1,N(j))
                !    W_N(2,1) = W_O(2,1)
                !    W_N(2,N(j)) = W_O(2,N(j))
                !else if (condition == 'P') then
                !    ! condition periodique
                !    W_N(:,1) = W_O(:,1) - (dt/dx)* (Flux(:,1) - Flux(:,(N(j)-1)))
                !   W_N(:,N(j)) = W_N(:,1)
                !else ! par defaut on prend des conditions de Neumann
                !    W_N(1,1) = W_N(1,2)
                !    W_N(1,N(j)) = W_N(1,N(j)-1)
                !    W_N(2,1) = W_N(2,2)
                !    W_N(2,N(j)) = W_N(2,N(j)-1)
                !end if

                ! conditions aux bords en amont
                do i = 1,2
                    !write(6,*) 'amont', i, cond(i,1)
                    if (cond(i,1) == -1.0_rp) then
                        W_N(i,1) = W_N(i,2)
                    else if (cond(i,1) == 0.0_rp) then
                        W_N(i,1) = W_O(i,1)
                    else
                        W_N(i,1) = cond(i,1)
                    end if
                end do
                ! conditions aux bords en aval
                do i = 1,2
                    !write(6,*) 'aval', i, cond(i,2)
                    if (cond(i,2) == -1.0_rp) then
                        W_N(i,N(j)) = W_N(i,N(j)-1)
                    else if (cond(i,2) == 0.0_rp) then
                        W_N(i,N(j)) = W_O(i,N(j))
                    else
                        W_N(i,1) = cond(i,2)
                    end if
                end do
                
                !mise a jour
                W_O(1:2,1:N(j)) = W_N(1:2,1:N(j))

                ! on calcule le nombre d'iterations 
                Nb_iter = Nb_iter + 1

                ! on etait en variables conservatives donc on passe en variables primitives pour le calcul de la CFL, ou pour la sauvegarde des donnees si derniere iteration
                call conserv_to_prim(W_O, N(j))
            end do

            write(6,*) 'Nombre d iterations: ', Nb_iter
            Err(2,j) = norme_L2(W_O(2,:),N(j))
            do i = 1,N(j)
                Err_h(i) = W_O(1,i) + Zi(i) - hR
            end do
            Err(1,j) = norme_L2(Err_h,N(j))
            write(6,*)
            Err(3,j) = dx

            deallocate(W_O, W_N, Flux, Err_u, Err_h, Zi, W_Om, W_Op)
        end do

        ! Enregistrement des erreurs
        write(6,*) 'Enregistrement dans le fichier erreurs.dat'
        call sauvegarde_conv('erreurs.dat', nb, Err)

        deallocate(Err, N)

    else ! si on ne veut pas faire de graphiques de convergence
        write(6,*) 'Calcul de la solution pour Ns cellules avec Ns = ', Ns
        write(6,*)
        ! calcul pas de maillage en espace dx
        dx = (x_fin - x_deb)/Ns

        ! allocation memoire des tableaux
        allocate(W_O(2,1:Ns), W_N(2,1:Ns), Flux(2,1:(Ns-1)), Err_u(Ns), Err_h(Ns), Zi(Ns))
        allocate(W_Om(2,1:Ns),W_Op(2,1:Ns))

        ! calcul topographie
        call topographie(Zi, Ns, dx, x_deb, topo)

        ! initialisation pour t = 0
        call initialisation_syst(W_O, Ns, x_deb, x_fin, uL, uR, hL, hR, Zi) ! en variables primitives

        ! boucle en temps
        date = 0._rp
        do while (date < T_fin)
            ! CFL 
            v = abs(lambda_1(W_O(:,1))) ! on fait le max des valeurs propres
            v = max(v,abs(lambda_2(W_O(:,1))))
            do i = 2,Ns
                v = max(v, abs(lambda_1(W_O(:,i))), abs(lambda_2(W_O(:,i))))
            end do
            dt = dx*CFL/v
            dt = min(dt, T_fin - date)
            date = date + dt

            ! on etait en variables primitives donc on passe en variables conservatives
            call prim_to_conserv(W_O, Ns)
            
            ! calcul des flux
            if (schema == 'LF') then
                call flux_LF_syst(Ns, Flux, W_O, dt, dx)
            else if (schema == 'RS') then
                call flux_RS_syst(Ns, Flux, W_O)
            else if (schema == 'HL') then
                call flux_HLL_syst(Ns, Flux, W_O, dx, dt)
            else if (schema == 'HY') then
                call flux_recons_hydro(Ns, Flux, W_O, Zi, dx, dt, W_Om, W_Op)
            end if

            ! update calcul de u_i^{n+1}
            Delta = dt/dx
            if (schema == 'HY') then
                do i = 2,(Ns-1)
                    W_N(1,i) = W_O(1,i) - Delta*(Flux(1,i) - Flux(1,i-1))
                    W_N(2,i) = W_O(2,i) - Delta*(Flux(2,i) - Flux(2,i-1)) + dt*terme_src_hy(dx,W_Op(1,i-1),W_Om(1,i))
                end do
            else
                do i = 2,(Ns-1)
                    W_N(1,i) = W_O(1,i) - Delta*(Flux(1,i) - Flux(1,i-1))
                    W_N(2,i) = W_O(2,i) - Delta*(Flux(2,i) - Flux(2,i-1)) + dt*terme_src(W_O(1,i),dx,Zi(i-1), Zi(i+1))
                end do
            end if

            ! Conditions aux limites
            !if (condition == 'D') then 
            !    ! Dirichlet
            !    W_N(1,1) = W_O(1,1)
            !    W_N(1,Ns) = W_O(1,Ns)
            !    W_N(2,1) = W_O(2,1)
            !    W_N(2,Ns) = W_O(2,Ns)
            !else if (condition == 'P') then
            !    ! condition periodique
            !    W_N(:,1) = W_O(:,1) - (dt/dx)* (Flux(:,1) - Flux(:,(Ns-1)))
            !    W_N(:,Ns) = W_N(:,1)
            !else ! par defaut on prend des conditions de Neumann
            !    W_N(1,1) = W_N(1,2)
            !    W_N(1,Ns) = W_N(1,Ns-1)
            !    W_N(2,1) = W_N(2,2)
            !    W_N(2,Ns) = W_N(2,Ns-1)
            !end if

            ! conditions aux bords en amont
            do i = 1,2
                !write(6,*) 'amont', i, cond(i,1)
                if (cond(i,1) == -1.0_rp) then
                    W_N(i,1) = W_N(i,2)
                else if (cond(i,1) == 0.0_rp) then
                    W_N(i,1) = W_O(i,1)
                else
                    W_N(i,1) = cond(i,1)
                end if
            end do
            ! conditions aux bords en aval
            do i = 1,2
                !write(6,*) 'aval', i, cond(i,2)
                if (cond(i,2) == -1.0_rp) then
                    W_N(i,Ns) = W_N(i,Ns-1)
                else if (cond(i,2) == 0.0_rp) then
                    W_N(i,Ns) = W_O(i,Ns)
                else
                    W_N(i,Ns) = cond(i,2)
                end if
            end do
            
            !mise a jour
            W_O(1:2,1:Ns) = W_N(1:2,1:Ns)

            ! on calcule le nombre d'iterations 
            Nb_iter = Nb_iter + 1

            ! on etait en variables conservatives donc on passe en variables primitives pour le calcul de la CFL, ou pour la sauvegarde des donnees si derniere iteration
            call conserv_to_prim(W_O, Ns)
        end do

        write(6,*) 'Nombre d iterations', Nb_iter
        write(6,*)
        write(6,*) 'Norme L^2 de u_i^n: ', norme_L2(W_O(2,:),Ns)
        do i = 1,Ns
            Err_h(i) = W_O(1,i) + Zi(i) - hR
        end do
        write(6,*) 'Norme L^2 de h_i^n+z_i-H: ', norme_L2(Err_h,Ns)
        write(6,*)
        
        ! on sauvegarde les resultats pour t = T_fin
        write(6,*) 'Enregistrement dans les fichiers solution_h.dat et solution_u.dat'
        write(6,*) 'Enregistrement de la topographie dans le fichier topo.dat'
        call sauvegarde_syst('solution_h.dat','solution_u.dat', W_O, Ns, x_deb, x_fin, Zi)

        deallocate(W_O, W_N, Flux, Err_u, Err_h, Zi, W_Om, W_Op)
    end if 

end program systeme_SW