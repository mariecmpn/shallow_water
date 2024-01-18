module initialisation_sauvegarde
    use numerics
    IMPLICIT NONE

    contains 

    ! Fonctions conditions initiales dans le cas d'un choc/contact (u_R < u_L)

    real(rp) function choc_u0(x)
        real(rp) :: x
        if (x<0.) then
            choc_u0 = 3._rp
        else 
            choc_u0 = 1._rp
        end if
    end function choc_u0

    real(rp) function choc_rho0(x)
        real(rp) :: x
        if (x<0.) then
            choc_rho0 = 1._rp
        else 
            choc_rho0 = 2._rp
        end if
    end function choc_rho0

    ! Fonctions conditions initiales dans le cas d'une detente/contact (u_R > u_L)

    real(rp) function det_u0(x)
        real(rp) :: x
        if (x<0.) then
            det_u0 = 1._rp
        else 
            det_u0 = 3._rp
        end if
    end function det_u0

    real(rp) function det_rho0(x)
        real(rp) :: x
        if (x<0.) then
            det_rho0 = 2._rp
        else 
            det_rho0 = 1._rp
        end if
    end function det_rho0

    ! Fonction condition initiale qui prend en compte la fonction demandee

    real(rp) function initial_u(x,fonc)
        ! fonction condition initiale pour u
        real(rp) :: x
        integer :: fonc
        if (fonc == 1) then
            initial_u = choc_u0(x)
        else 
            initial_u = det_u0(x)
        end if
    end function initial_u


    real(rp) function initial_rho(x, fonc)
        ! fonction condition initiale pour rho
        real(rp) :: x
        integer :: fonc
        if (fonc == 1) then
            initial_rho = choc_rho0(x)
        else 
            initial_rho = det_rho0(x)
        end if
    end function initial_rho

    ! Fonction exacte

    function sol_ex_syst(x,t,fonc,v_max,rho_max)
        ! solutions exactes du systeme en fonction de la fonction choisie
        real(rp), dimension(2) :: sol_ex_syst ! vecteur W = (rho, u) solution exacte
        real(rp) :: v_max,rho_max ! parametres de la fonction p
        integer :: fonc ! fonction choisie
        real(rp) :: x,t ! coordonnees x en espace et t en temps
        real(rp) :: sigma, u, rho, ul, ur, rhol, rhor, rhoe

        if (fonc == 1) then ! fonction choc/contact
            ul = initial_u(-1._rp,1) ! on definit u_L, u_R, rho_L et rho_R
            ur = initial_u(1._rp,1)
            rhol = initial_rho(-1._rp,1)
            rhor = initial_rho(1._rp,1)
            rho = sqrt(rhol**2 + rho_max/v_max*(ul-ur))
            u = ur
            sigma = (rho*u - rhol*ul)/ (rho - rhol)
            if (x<sigma*t) then ! avant le choc
                sol_ex_syst(1) = rhol
                sol_ex_syst(2) = ul
            else if (x>=sigma*t .AND. x<ur*t) then ! entre choc et contact
                sol_ex_syst(1) = rho
                sol_ex_syst(2) = ur
            else ! apres contact
                sol_ex_syst(1) = rhor
                sol_ex_syst(2) = ur
            end if 
        else if (fonc == 2) then ! fonction detente/contact
            ul = initial_u(-1._rp,2) ! on definit u_L, u_R, rho_L et rho_R
            ur = initial_u(1._rp,2)
            rhol = initial_rho(-1._rp,2)
            rhor = initial_rho(1._rp,2)
            rhoe = sqrt(rhol**2 + (rho_max/v_max)*(ul-ur)) ! rho^*
            !write(6,*) (ul-rhol*p_prime(rhol,v_max,rho_max))*t
            !write(6,*) (ur-rhor*p_prime(rhor,v_max,rho_max))*t
            !write(6,*) x/t
            !write(6,*)
            if (x<(ul-rhol*p_prime(rhol,v_max,rho_max))*t) then
                sol_ex_syst(1) = rhol
                sol_ex_syst(2) = ul
            else if (x>=(ul-rhol*p_prime(rhol,v_max,rho_max))*t .AND. x<(ur-rhor*p_prime(rhor,v_max,rho_max))*t) then
                sol_ex_syst(1) = sqrt(rho_max*(ul - x/t)/(3.*v_max)+((rhol)**2)/3.)
                sol_ex_syst(2) = ul-p(sol_ex_syst(1),v_max,rho_max)+p(rhol,v_max,rho_max)
            else if (x>=(ur-rhor*p_prime(rhor,v_max,rho_max))*t .AND. x<ur*t) then
                sol_ex_syst(1) = sqrt(rhol**2 + (rho_max/v_max)*(ul-ur))
                sol_ex_syst(2) = ur
            else
                sol_ex_syst(1) = rhor
                sol_ex_syst(2) = ur
            end if
        end if 
    end function sol_ex_syst

    ! Lecture des donnees, initialisation et sauvegarde

    subroutine lecture_donnees_syst(file_name, x_deb, x_fin, Ns, CFL, T_fin, condition, schema, v_max, rho_max, fonc)
        ! Subroutine pour recuperer les donnees du fichier file_name
        IMPLICIT NONE
        character(len = *), intent(in) :: file_name ! nom du fichier a ouvrir
        integer, intent(inout) :: Ns ! nombre de cellules
        real(rp), intent(inout) :: x_deb, x_fin ! debut et fin des x
        real(rp), intent(inout) :: CFL ! condition CFL
        real(rp), intent(inout) :: T_fin ! temps final
        character(len = 1), intent(inout) :: condition ! condition aux bords
        integer, intent(inout) :: fonc ! fonction initiale utilisee
        character(len = 2), intent(inout) :: schema ! schema utilise
        real(rp), intent(inout) :: v_max, rho_max ! vitesse max et densite max pour p

        integer :: my_unit

        open(newunit = my_unit, file = file_name, action = 'read', form = 'formatted', status = 'old')
        
        read(my_unit, *) x_deb, x_fin
        read(my_unit, *) Ns
        read(my_unit, *) CFL
        read(my_unit, *) T_fin
        read(my_unit, *) condition
        read(my_unit, *) schema
        read(my_unit, *) rho_max
        read(my_unit, *) v_max
        read(my_unit, *) fonc

        !if (condition_lim == 'D') then
        !    condition = 0
        !else if (condition_lim == 'N') then
        !    condition = 1
        !else
        !    condition = -1
        !end if

        if (schema == 'LF') then ! Lax-Friedrichs
            write(6,*) "Schema utilise: Lax_Friedrichs"
        else if (schema == 'RS') then ! Rusanov
            write(6,*) "Schema utilise: Rusanov"
        else if (schema == 'HL') then
            write(6,*) "Schema utilise: HLL"
        !else if (schema_use == 'GD') then ! Godunov
        !    schema = 2
        !else if (schema_use == 'LW') then ! Lax-Wendroff
        !    schema = 3
        !else
        !    schema = 0
        end if

        close(my_unit)
    end subroutine lecture_donnees_syst


    subroutine initialisation_syst(W_O, Ns, x_deb, x_fin, fonc)
        ! suboutine pour initialiser le probleme
        IMPLICIT NONE
        integer, intent(in) :: Ns, fonc ! nmbre de cellules, fonction utilisee
        real(rp), dimension(2,1:Ns), intent(inout) :: W_O ! tableau qu'on initialise
        real(rp), intent(in) :: x_deb, x_fin ! debut et fin des x
        real(rp) :: x
        integer :: i
        real(rp) :: deltax

        deltax = (x_fin-x_deb)/Ns
        do i = 1,Ns
            x = x_deb + i*deltax
            W_O(1,i) = initial_rho(x, fonc)
            W_O(2,i) = initial_u(x, fonc)
        end do
    end subroutine initialisation_syst


    subroutine sauvegarde_syst(file_name_rho, file_name_u, W_O, Ns, x_deb, x_fin)
        ! subroutine pour sauvegarde les solutions du probleme
        IMPLICIT NONE
        character(len = *), intent(in) :: file_name_u, file_name_rho ! noms des fichiers de sortie
        integer, intent(in) :: Ns ! nombre de cellules
        real(rp), dimension(2,1:Ns), intent(in) :: W_O ! tableau a enregistrer
        real(rp), intent(in) :: x_deb, x_fin ! debut et fin des x
        real(rp) :: x ! pour calculer x_i
        integer :: i ! pour boucle do
        integer :: my_unit_1 = 60, my_unit_2 = 70

        open(my_unit_1, file = file_name_rho, action = 'write', form = 'formatted', status = 'unknown')
        open(my_unit_2, file = file_name_u, action = 'write', form = 'formatted', status = 'unknown')

        do i = 1,Ns
            x = x_deb + i*(x_fin-x_deb)/Ns
            write(my_unit_1, *) x, W_O(1,i)
            write(my_unit_2, *) x, W_O(2,i)
        end do

        close(my_unit_1)
        close(my_unit_2)
    end subroutine sauvegarde_syst

    ! Fonction p, p' et calcule de la vitesse (en prenant en compte rho = 0)

    real(rp) function p(rho, v_max, rho_max)
        ! fonction p
        real(rp) :: rho, v_max, rho_max
        p = v_max*(rho**2/rho_max)
    end function p


    real(rp) function p_prime(rho, v_max, rho_max)
        ! fonction derivee de p
        real(rp) :: rho, v_max, rho_max
        p_prime = 2._rp*v_max*(rho/rho_max)
    end function p_prime


    real(rp) function vitesse(W, v_max, rho_max)
        ! fonction pour la vitesse
        ! au cas ou rho = 0
        real(rp) :: eps = 10.E-14
        real(rp), dimension(2), intent(in) :: W
        real(rp) :: rho, u
        real(rp) :: v_max, rho_max

        rho = W(1)
        if (rho > eps) then 
            u = (W(2) / rho) - p(rho, v_max, rho_max)
        else
            u = 0
        end if
        vitesse = u 
    end function vitesse

    ! Pour passer des variables conservatives aux variables primitives et inversement

    subroutine conserv_to_prim(W_O, Ns, v_max, rho_max)
        ! Subroutine pour passer des variables conservatives aux variables primitives
        ! Q = rho*u + rho*p(rho) devient u
        integer, intent(in) :: Ns
        real(rp), dimension(2,1:Ns), intent(inout) :: W_O
        real(rp), intent(in) :: v_max, rho_max
        integer :: i
        do i = 1,Ns
            W_O(2,i) = vitesse(W_O(1:2,i), v_max, rho_max)
        end do
    end subroutine conserv_to_prim


    subroutine prim_to_conserv(W_O, Ns, v_max, rho_max)
        ! Subroutine pour passer des variables primitives aux variables conservatives
        ! u devient Q = rho*u + rho*p(rho)
        integer, intent(in) :: Ns
        real(rp), dimension(2,1:Ns), intent(inout) :: W_O
        real(rp), intent(in) :: v_max, rho_max
        integer :: i
        do i = 1,Ns
            W_O(2,i) = W_O(1,i)*W_O(2,i) + W_O(1,i)*p(W_O(1,i), v_max, rho_max)
        end do
    end subroutine prim_to_conserv


    ! Algebre du systeme (matrice A, fonction F...)

    real(rp) function lambda_1(W, v_max, rho_max)
        ! fonction pour la valeur propre 1 du systeme
        ! En variables primitives
        real(rp), dimension(2) :: W
        real(rp) :: v_max, rho_max
        lambda_1 = W(2) - W(1)*p_prime(W(1), v_max, rho_max)
    end function lambda_1


    real(rp) function lambda_2(W)
        ! fonction pour la valeur propre 1 du systeme
        ! En variables primitives
        real(rp), dimension(2) :: W 
        lambda_2 = W(2)
    end function lambda_2


    function mat_A(u, rho, v_max, rho_max)
        ! fonction qui calcule la matrice A
        real(rp), intent(in) :: u, rho, v_max, rho_max
        real(rp), dimension(2,2) :: mat_A
        mat_A(1,1) = u
        mat_A(1,2) = rho
        mat_A(2,1) = 0._rp
        mat_A(2,2) = u - rho*p_prime(rho, v_max, rho_max)
    end function mat_A


    subroutine F(F_ex, W, v_max, rho_max)
        ! subroutine pour calculer la fonction F flux exact
        ! W en variables conservatives
        real(rp), dimension(2), intent(in) :: W
        real(rp), dimension(2), intent(inout) :: F_ex
        real(rp) ,intent(in) ::  v_max, rho_max
        real(rp) :: u
        u = vitesse(W(1:2), v_max, rho_max) ! on calcule u car on est en variables conservatives
        F_ex(:) = u * W(:) ! on multiplie par u car u est 
    end subroutine F

    real(rp) function norme_L2(U, Ns)
        integer :: Ns
        real(rp), dimension(Ns) :: U
        integer :: i
        norme_L2 = 0._rp
        do i = 1,Ns
            norme_L2 = norme_L2+U(i)**2
        end do
        norme_L2 = sqrt(norme_L2)
    end function norme_L2


end module initialisation_sauvegarde