module initialisation_sauvegarde
    use numerics
    IMPLICIT NONE

    contains 

    ! Fonction condition initiale qui prend en compte la fonction demandee

    real(rp) function initial_u(x, uL, uR)
        ! fonction condition initiale pour u
        real(rp) :: x, uL, uR
        if (x < 0.) then
            initial_u = uL
        else
            initial_u = uR
        end if
    end function initial_u


    real(rp) function initial_h(x, hL, hR)
        ! fonction condition initiale pour rho
        real(rp) :: x, hL, hR
        if (x < 0.) then
            initial_h = hL
        else
            initial_h = hR
        end if
    end function initial_h

    ! Lecture des donnees, initialisation et sauvegarde

    subroutine lecture_donnees_syst(file_name,x_deb,x_fin,Ns,CFL,T_fin,condition,schema,uL,uR,hL,hR,topo)
        ! Subroutine pour recuperer les donnees du fichier file_name
        IMPLICIT NONE
        character(len = *), intent(in) :: file_name ! nom du fichier a ouvrir
        integer, intent(inout) :: Ns ! nombre de cellules
        real(rp), intent(inout) :: x_deb, x_fin ! debut et fin des x
        real(rp), intent(inout) :: CFL ! condition CFL
        real(rp), intent(inout) :: T_fin ! temps final
        character(len = 1), intent(inout) :: condition ! condition aux bords
        real(rp), intent(inout) :: uL, uR, hL, hR ! conditions initiales
        character(len = 2), intent(inout) :: schema ! schema utilise
        integer, intent(inout) :: topo ! si on veut une topographie ou non, et laquelle si oui

        integer :: my_unit

        open(newunit = my_unit, file = file_name, action = 'read', form = 'formatted', status = 'old')
        
        read(my_unit, *) x_deb, x_fin
        read(my_unit, *) Ns
        read(my_unit, *) CFL
        read(my_unit, *) T_fin
        read(my_unit, *) condition
        read(my_unit, *) schema
        read(my_unit, *) hL, hR
        read(my_unit, *) uL, uR
        read(my_unit, *) topo

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

        write(6,*)
        if (topo == 0) then
            write(6,*) 'Pas de topographie'
        else if (topo == 1) then
            write(6,*) 'Topographie choisie: Z(x) = (0.2-0.05(x-10)^2)_+'
        end if

        close(my_unit)
    end subroutine lecture_donnees_syst


    subroutine initialisation_syst(W_O, Ns, x_deb, x_fin, uL, uR, hL, hR, Zi)
        ! suboutine pour initialiser le probleme
        IMPLICIT NONE
        integer, intent(in) :: Ns ! nmbre de cellules, fonction utilisee
        real(rp), dimension(2,1:Ns), intent(inout) :: W_O ! tableau qu'on initialise
        real(rp), intent(in) :: x_deb, x_fin ! debut et fin des x
        real(rp), intent(in) :: uL, uR, hL, hR
        real(rp), dimension(Ns), intent(in) :: Zi
        real(rp) :: x
        integer :: i
        real(rp) :: deltax

        deltax = (x_fin-x_deb)/Ns
        do i = 1,Ns
            x = x_deb + i*deltax
            W_O(1,i) = initial_h(x, hL, hR) - Zi(i)
            W_O(2,i) = initial_u(x, uL, uR)
        end do
    end subroutine initialisation_syst


    subroutine sauvegarde_syst(file_name_h, file_name_u, W_O, Ns, x_deb, x_fin, Zi, topo)
        ! subroutine pour sauvegarde les solutions du probleme
        IMPLICIT NONE
        character(len = *), intent(in) :: file_name_u, file_name_h ! noms des fichiers de sortie
        integer, intent(in) :: Ns ! nombre de cellules
        real(rp), dimension(2,1:Ns), intent(in) :: W_O ! tableau a enregistrer
        real(rp), dimension(Ns), intent(in) :: Zi ! topographie a enregistrer si non nulle
        real(rp), intent(in) :: x_deb, x_fin ! debut et fin des x
        integer, intent(in) :: topo
        real(rp) :: x ! pour calculer x_i
        integer :: i ! pour boucle do
        integer :: my_unit_1 = 60, my_unit_2 = 70, my_unit_3 = 80

        open(my_unit_1, file = file_name_h, action = 'write', form = 'formatted', status = 'unknown')
        open(my_unit_2, file = file_name_u, action = 'write', form = 'formatted', status = 'unknown')
        if (topo == 1) open(my_unit_3, file = 'topo.dat', action = 'write', form = 'formatted', status = 'unknown')

        do i = 1,Ns
            x = x_deb + i*(x_fin-x_deb)/Ns
            write(my_unit_1, *) x, W_O(1,i)+Zi(i)
            write(my_unit_2, *) x, W_O(2,i)
            if (topo == 1) write(my_unit_3, *) x, Zi(i)
        end do

        close(my_unit_1)
        close(my_unit_2)
        close(my_unit_3)
    end subroutine sauvegarde_syst

    subroutine sauvegarde_topo(file_name, Zi, Ns, x_deb, x_fin)
        character(len = *), intent(in) :: file_name
        integer, intent(in) :: Ns ! nombre de cellules
        real(rp), dimension(Ns), intent(in) :: Zi ! topographie a enregistrer
        real(rp), intent(in) :: x_deb, x_fin ! debut et fin des x
        real(rp) :: x ! pour calculer x_i
        integer :: i ! pour boucle do
        integer :: my_unit

        open(newunit = my_unit, file = file_name, action = 'write', form = 'formatted', status = 'unknown')

        do i = 1,Ns
            x = x_deb + i*(x_fin-x_deb)/Ns
            write(my_unit, *) x, Zi(i)
        end do

        close(my_unit)
    end subroutine sauvegarde_topo

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

    real(rp) function Z(x) ! fonction pour la topographie
        real(rp) :: x, y
        y = 0.2_rp - 0.05_rp*(x-10.0_rp)**2
        if (y > 0.) then 
            Z = y
        else
            Z = 0._rp
        end if
    end function Z

    real(rp) function terme_src(h, dx, zi, zj)
        real(rp) :: h, dx
        real(rp) :: zi,zj

        terme_src = -g*h*(zj-zi)/(2._rp*dx)
    end function terme_src

    subroutine topographie(Zi, Ns, dx, x_deb, topo)
        real(rp), dimension(Ns), intent(out) :: Zi
        integer, intent(in) :: Ns, topo
        real(rp), intent(in) :: dx
        real(rp), intent(in) :: x_deb
        integer :: i

        Zi(:) = 0._rp
        if (topo == 1) then
            do i = 1,Ns
                Zi(i) = Z(x_deb+i*dx)
            end do
        end if

    end subroutine topographie


end module initialisation_sauvegarde