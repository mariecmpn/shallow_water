program systeme_trafic
    use numerics
    use initialisation_sauvegarde
    use schemasSH
    IMPLICIT NONE
    real(rp) :: x_deb, x_fin
    integer :: Ns
    real(rp) :: CFL, T_fin, date, dt, dx
    real(rp) :: v, Delta, x
    integer :: i
    real(rp), dimension(:,:), allocatable :: W_O
    real(rp), dimension(:,:), allocatable :: W_N
    real(rp), dimension(:,:), allocatable :: Flux
    real(rp), dimension(:,:), allocatable :: W_ex
    real(rp), dimension(:), allocatable :: Err_u, Err_rho
    character(len = 1) :: condition
    character(len = 2) :: schema
    integer :: fonc
    integer :: Nb_iter = 0

    write(6,*) '------------------------------------------'
    write(6,*) '----------- Resolution syteme ------------'
    write(6,*) '------------- SHALLOW-WATER --------------'
    write(6,*) '------------------------------------------'

    write(6,*)

    ! lecture des donnees du fichier donnees.dat
    call lecture_donnees_syst('init.dat', x_deb, x_fin, Ns, CFL, T_fin, condition, schema, fonc)

    write(6,*)

    dx = (x_fin - x_deb)/Ns

    ! allocation memoire des tableaux
    allocate(W_O(2,1:Ns), W_N(2,1:Ns), Flux(2,1:(Ns-1)), W_ex(2,1:Ns), Err_u(Ns), Err_rho(Ns))

    ! initialisation pour t = 0
    call initialisation_syst(W_O, Ns, x_deb, x_fin, fonc) ! en variables primitives

    ! boucle en temps
    date = 0._rp
    do while (date < T_fin)
        ! CFL (raccourci u_i^n + rho_i^n*p'(rho_i^n) pour max des valeurs propres)
        v = abs(lambda_1(W_O(:,1)))
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
            W_N(2,i) = W_O(2,i) - Delta*(Flux(2,i) - Flux(2,i-1))
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
    write(6,*) 'Nombre de cellules: ', Ns
    !write(6,*) 'Erreurs L2 entre solution approchee et solution exacte: ' 
    !write(6,*) 'Pour rho: ', norme_L2(Err_rho, Ns)
    !write(6,*) 'Pour u: ', norme_L2(Err_u, Ns)
    !write(6,*) 'Erreurs L2 des solutions approchees: ' 
    !write(6,*) 'Pour rho: ', norme_L2(W_O(1,:), Ns)
    !write(6,*) 'Pour u: ', norme_L2(W_O(2,:), Ns)
    write(6,*)
    
    ! on sauvegarde les resultats pour t = T_fin
    call sauvegarde_syst('solution_h.dat','solution_u.dat', W_O, Ns, x_deb, x_fin)

    ! on sauvegarde la fonction exacte pour t = T_fin
    !call sauvegarde_syst('solution_rho_ex.dat', 'solution_u_ex.dat', W_ex, Ns, x_deb, x_fin)

    deallocate(W_O, W_N, Flux, Err_u, Err_rho, W_ex)

end program systeme_trafic