program systeme_SW
    use numerics
    use initialisation_sauvegarde
    use schemasSW
    IMPLICIT NONE
    real(rp) :: x_deb, x_fin
    integer :: Ns
    real(rp) :: CFL, T_fin, date, dt, dx
    real(rp) :: v, Delta, x
    integer :: i
    real(rp) :: uL, uR, hL, hR
    real(rp), dimension(:,:), allocatable :: W_O
    real(rp), dimension(:,:), allocatable :: W_N
    real(rp), dimension(:,:), allocatable :: Flux
    real(rp), dimension(:), allocatable :: Zi
    real(rp), dimension(:), allocatable :: Err_u, Err_h
    character(len = 1) :: condition
    character(len = 2) :: schema
    integer :: Nb_iter = 0
    integer :: topo

    write(6,*) '------------------------------------------'
    write(6,*) '----------- Resolution syteme ------------'
    write(6,*) '------------- SHALLOW-WATER --------------'
    write(6,*) '------------------------------------------'

    write(6,*)

    ! lecture des donnees du fichier donnees.dat
    call lecture_donnees_syst('init.dat', x_deb, x_fin, Ns, CFL, T_fin, condition, schema, uL, uR, hL, hR, topo)

    write(6,*)

    dx = (x_fin - x_deb)/Ns

    ! allocation memoire des tableaux
    allocate(W_O(2,1:Ns), W_N(2,1:Ns), Flux(2,1:(Ns-1)), Err_u(Ns), Err_h(Ns), Zi(Ns))

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
            call flux_HLL_syst(Ns, Flux, W_O)
        !else if (schema == 1) then
        !    call flux_MR(Ns, Flux, W_O)
        !else if (schema == 2) then
        !    call flux_GD(Ns, Flux, W_O)
        !else if (schema == 3) then
        !   call flux_LW(Ns, Flux, W_O, dt, dx)
        end if

        ! update calcul de u_i^{n+1}
        Delta = dt/dx
        do i = 2,(Ns-1)
            W_N(1,i) = W_O(1,i) - Delta*(Flux(1,i) - Flux(1,i-1))
            W_N(2,i) = W_O(2,i) - Delta*(Flux(2,i) - Flux(2,i-1)) + dt*terme_src(W_O(1,i),dx,Zi(i-1), Zi(i+1))
        end do

        ! Conditions aux limites
        if (condition == 'D') then 
            ! Dirichlet
            W_N(1,1) = W_O(1,1)
            W_N(1,Ns) = W_O(1,Ns)
            W_N(2,1) = W_O(2,1)
            W_N(2,Ns) = W_O(2,Ns)
        else if (condition == 'P') then
            ! condition periodique
            W_N(:,1) = W_O(:,1) - (dt/dx)* (Flux(:,1) - Flux(:,(Ns-1)))
            W_N(:,Ns) = W_N(:,1)
        else ! par defaut on prend des conditions de Neumann
            W_N(1,1) = W_N(1,2)
            W_N(1,Ns) = W_N(1,Ns-1)
            W_N(2,1) = W_N(2,2)
            W_N(2,Ns) = W_N(2,Ns-1)
        end if
        
        !mise a jour
        W_O(1:2,1:Ns) = W_N(1:2,1:Ns)

        ! on calcule le nombre d'iterations 
        Nb_iter = Nb_iter + 1

        ! on etait en variables conservatives donc on passe en variables primitives pour le calcul de la CFL, ou pour la sauvegarde des donnees si derniere iteration
        call conserv_to_prim(W_O, Ns)
    end do

    write(6,*) 'Nombre d iterations', Nb_iter
    write(6,*)
    write(6,*) 'Norme L^2 de u_i^n: ', norme_L2(W_O(1,:),Ns)
    do i = 1,Ns
        Err_h(i) = W_O(1,i) + Zi(i) - hR
    end do
    write(6,*) 'Norme L^2 de h_i^n+z_i-H: ', norme_L2(Err_h,Ns)
    write(6,*)
    
    ! on sauvegarde les resultats pour t = T_fin
    write(6,*) 'Enregistrement dans les fichiers solution_h.dat et solution_u.dat'
    if (topo == 1) write(6,*) 'Enregistrement de la topographie dans le fichier topo.dat'
    call sauvegarde_syst('solution_h.dat','solution_u.dat', W_O, Ns, x_deb, x_fin, Zi, topo)


    deallocate(W_O, W_N, Flux, Err_u, Err_h)

end program systeme_SW