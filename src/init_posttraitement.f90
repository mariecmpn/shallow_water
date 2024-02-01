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

    real(rp) function choc_h0(x)
        real(rp) :: x
        if (x<0.) then
            choc_h0 = 1._rp
        else 
            choc_h0 = 2._rp
        end if
    end function choc_h0

    ! Fonctions conditions initiales dans le cas d'une detente/contact (u_R > u_L)

    real(rp) function det_u0(x)
        real(rp) :: x
        if (x<0.) then
            det_u0 = 1._rp
        else 
            det_u0 = 3._rp
        end if
    end function det_u0

    real(rp) function det_h0(x)
        real(rp) :: x
        if (x<0.) then
            det_h0 = 2._rp
        else 
            det_h0 = 1._rp
        end if
    end function det_h0

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


    real(rp) function initial_h(x, fonc)
        ! fonction condition initiale pour rho
        real(rp) :: x
        integer :: fonc
        if (fonc == 1) then
            initial_h = choc_h0(x)
        else 
            initial_h = det_h0(x)
        end if
    end function initial_h

    ! Lecture des donnees, initialisation et sauvegarde

    subroutine lecture_donnees_syst(file_name, x_deb, x_fin, Ns, CFL, T_fin, condition, schema, fonc)
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

        integer :: my_unit

        open(newunit = my_unit, file = file_name, action = 'read', form = 'formatted', status = 'old')
        
        read(my_unit, *) x_deb, x_fin
        read(my_unit, *) Ns
        read(my_unit, *) CFL
        read(my_unit, *) T_fin
        read(my_unit, *) condition
        read(my_unit, *) schema
        read(my_unit, *) fonc

        write(6,*) 'Calcul entre x_deb = ', x_deb, ' et x_fin = ', x_fin
        write(6,*) 'Nombre de points de maillage: ', Ns
        write(6,*) 'Temps final: ', T_fin

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
            W_O(1,i) = initial_h(x, fonc)
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

    ! Pour passer des variables conservatives aux variables primitives et inversement

    subroutine conserv_to_prim(W_O, Ns)
        ! Subroutine pour passer des variables conservatives aux variables primitives
        ! Q = hu devient u
        integer, intent(in) :: Ns
        real(rp), dimension(2,1:Ns), intent(inout) :: W_O
        integer :: i
        do i = 1,Ns
            W_O(2,i) = W_O(2,i)/W_O(1,i)
        end do
    end subroutine conserv_to_prim


    subroutine prim_to_conserv(W_O, Ns)
        ! Subroutine pour passer des variables primitives aux variables conservatives
        ! u devient Q = hu
        integer, intent(in) :: Ns
        real(rp), dimension(2,1:Ns), intent(inout) :: W_O
        integer :: i
        do i = 1,Ns
            W_O(2,i) = W_O(1,i)*W_O(2,i)
        end do
    end subroutine prim_to_conserv

    ! Algebre du systeme (matrice A, fonction F...)

    real(rp) function lambda_1(W)
        ! fonction pour la valeur propre 1 du systeme
        ! En variables primitives
        real(rp), dimension(2) :: W
        lambda_1 = W(2) - sqrt(g*W(1))
    end function lambda_1

    real(rp) function lambda_2(W)
        ! fonction pour la valeur propre 1 du systeme
        ! En variables primitives
        real(rp), dimension(2) :: W 
        lambda_2 = W(2) + sqrt(g*W(1))
    end function lambda_2

    real(rp) function lambda_1_cons(W)
        ! fonction pour la valeur propre 1 du systeme
        ! En variables conservatives
        real(rp), dimension(2) :: W
        lambda_1_cons = W(2)/W(1) - sqrt(g*W(1))
    end function lambda_1_cons

    real(rp) function lambda_2_cons(W)
        ! fonction pour la valeur propre 1 du systeme
        ! En variables conservatives
        real(rp), dimension(2) :: W
        lambda_2_cons = W(2)/W(1) + sqrt(g*W(1))
    end function lambda_2_cons

    function A(W)
        ! matrice A en variables primitives
        real(rp), dimension(2,2) :: A
        real(rp), dimension(2) :: W
        A(1,1) = W(2)
        A(1,2) = W(1)
        A(2,1) = g
        A(2,2) = W(2)
    end function A

    function grad_F(W)
        ! matrice gradient de F en variables conservatives
        real(rp), dimension(2,2) :: grad_F
        real(rp), dimension(2) :: W
        grad_F(1,1) = 0._rp
        grad_F(1,2) = 1._rp
        grad_F(2,1) = -(W(2)/W(1))**2+g*W(1)
        grad_F(2,2) = 2.*(W(2)/W(1))
    end function grad_F

    function F(W)
        ! fonction pour calculer la fonction F flux exact
        ! W en variables conservatives
        real(rp), dimension(2) :: W
        real(rp), dimension(2) :: F
        F(1) = W(2)
        F(2) = (W(2)**2)/W(1) + 0.5_rp*g*W(1)**2
    end function F

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