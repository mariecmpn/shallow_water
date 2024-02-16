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

    subroutine lecture_donnees_syst(file_name,x_deb,x_fin,Ns,CFL,T_fin,condition,schema,uL,uR,hL,hR,topo,conv,Ns1,nb)
        ! Subroutine pour recuperer les donnees du fichier file_name
        IMPLICIT NONE
        character(len = *), intent(in) :: file_name ! nom du fichier a ouvrir
        integer, intent(inout) :: Ns ! nombre de cellules (si graphes de convergence => nombre de cellules min)
        real(rp), intent(inout) :: x_deb, x_fin ! debut et fin des x
        real(rp), intent(inout) :: CFL ! condition CFL
        real(rp), intent(inout) :: T_fin ! temps final
        character(len = 1), intent(inout) :: condition ! condition aux bords
        real(rp), intent(inout) :: uL, uR, hL, hR ! conditions initiales
        character(len = 2), intent(inout) :: schema ! schema utilise
        character(len = 1), intent(inout) :: conv ! si on veut faire ou non des graphiques de convergence
        integer, intent(inout) :: Ns1, nb ! si 
        integer, intent(inout) :: topo ! si on veut une topographie ou non, et laquelle si oui

        integer :: my_unit

        open(newunit = my_unit, file = file_name, action = 'read', form = 'formatted', status = 'old')
        
        read(my_unit, *) x_deb, x_fin
        read(my_unit, *) conv
        if (conv == 'Y') then ! si on veut faire de graphique de convergence
            read(my_unit, *) Ns, Ns1, nb
        else
            read(my_unit, *) Ns
        end if
        read(my_unit, *) CFL
        read(my_unit, *) T_fin
        read(my_unit, *) condition
        read(my_unit, *) schema
        read(my_unit, *) hL, hR
        read(my_unit, *) uL, uR
        read(my_unit, *) topo

        write(6,*) 'Calcul entre x_deb = ', x_deb, ' et x_fin = ', x_fin
        write(6,*) 'Temps final: ', T_fin

        if (schema == 'LF') then ! Lax-Friedrichs
            write(6,*) "Schema utilise: Lax_Friedrichs"
        else if (schema == 'RS') then ! Rusanov
            write(6,*) "Schema utilise: Rusanov"
        else if (schema == 'HL') then ! HLL
            write(6,*) "Schema utilise: HLL"
        else if (schema == 'HY') then ! reconstruction hydrostatique
            write(6,*) "Schema utilise: Reconstruction hydrostatique (avec schema HLL)"
        end if

        write(6,*)
        if (topo == 2) then
            write(6,*) 'Topographie choisie: Z(x) = (1.8-0.5(x-10)^2)_+'
        else if (topo == 1) then
            write(6,*) 'Topographie choisie: Z(x) = (0.2-0.05(x-10)^2)_+'
        else 
            write(6,*) 'Pas de topographie'
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
        if (topo == 1 .OR. topo == 2) open(my_unit_3, file = 'topo.dat', action = 'write', form = 'formatted', status = 'unknown')

        do i = 1,Ns
            x = x_deb + i*(x_fin-x_deb)/Ns
            write(my_unit_1, *) x, W_O(1,i)+Zi(i)
            write(my_unit_2, *) x, W_O(2,i)
            if (topo == 1 .OR. topo == 2) write(my_unit_3, *) x, Zi(i)
        end do

        close(my_unit_1)
        close(my_unit_2)
        close(my_unit_3)
    end subroutine sauvegarde_syst

    subroutine sauvegarde_conv(file_name, nb, Err)
        character(len = *), intent(in) :: file_name
        integer, intent(in) :: nb ! nombre de cellules
        !integer, dimension(nb), intent(in) :: N
        real(rp), dimension(3,nb), intent(in) :: Err
        integer :: i ! pour boucle do
        integer :: my_unit

        open(newunit = my_unit, file = file_name, action = 'write', form = 'formatted', status = 'unknown')

        do i = 1,nb
            write(my_unit, *) Err(3,i), Err(1,i), Err(2,i)
        end do

        close(my_unit)
    end subroutine sauvegarde_conv

    ! Pour passer des variables conservatives aux variables primitives et inversement

    subroutine conserv_to_prim(W_O, Ns)
        ! Subroutine pour passer des variables conservatives aux variables primitives
        ! Q = hu devient u
        integer, intent(in) :: Ns
        real(rp), dimension(2,1:Ns), intent(inout) :: W_O
        integer :: i
        do i = 1,Ns
            W_O(2,i) = vitesse(W_O(:,i))
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

    real(rp) function vitesse(W)
        real(rp), dimension(2) :: W
        if (W(1) > 1.E-6) then 
            vitesse = W(2)/W(1)
        else 
            vitesse = 0._rp
        end if

    end function vitesse

    real(rp) function Z(x,topo) ! fonction pour la topographie
        real(rp) :: x
        integer :: topo
        if (topo == 1) then 
            Z = max(0._rp, 0.2_rp - 0.05_rp*(x-10.0_rp)**2)
        else if (topo == 2) then
            Z = max(0._rp, 1.8_rp - 0.5_rp*(x-10.0_rp)**2)
        else 
            Z = 0._rp
        end if
        
    end function Z

    real(rp) function terme_src(h, dx, zi, zj)
        ! fonction pour le calcul du terme source
        real(rp) :: h, dx
        real(rp) :: zi,zj

        terme_src = -g*h*(zj-zi)/(2._rp*dx)
    end function terme_src

    real(rp) function terme_src_hy(dx, hi, hj)
        ! fonction pour le calcul du terme source pour la reconstruction hydrostatique
        real(rp) :: dx
        real(rp) :: hi,hj

        terme_src_hy = g*(hj**2-hi**2)/(2._rp*dx)
    end function terme_src_hy

    subroutine topographie(Zi, Ns, dx, x_deb, topo)
        ! subroutine pour le calcul du tableau de la topographie
        real(rp), dimension(Ns), intent(out) :: Zi
        integer, intent(in) :: Ns, topo
        real(rp), intent(in) :: dx
        real(rp), intent(in) :: x_deb
        integer :: i

        Zi(:) = 0._rp
        do i = 1,Ns
            Zi(i) = Z(x_deb+i*dx, topo)
        end do
    end subroutine topographie


end module initialisation_sauvegarde