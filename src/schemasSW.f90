module schemasSW
    use numerics
    use initialisation_sauvegarde
    IMPLICIT NONE

    contains

    subroutine flux_LF_syst(Ns, Flux, W_O, dt, dx)
    ! FLUX POUR SCHEMA DE LAX-FRIEDRICHS
        IMPLICIT NONE
        real(rp), intent(in) :: dt, dx
        integer, intent(in) :: Ns
        real(rp), dimension(2,Ns), intent(inout) :: Flux
        real(rp), dimension(2,Ns), intent(in) :: W_O
        !real(rp), dimension(2) :: F_ex1, F_ex
        real(rp) :: Delta
        integer :: i

        ! calcul flux pour Lax-Friedrichs
        Delta = (dx / dt) * 0.5_rp
        do i = 1,(Ns-1)
            Flux(:,i) = 0.5_rp*(F(W_O(:,i)) + F(W_O(:,i+1))) - Delta * (W_O(:,i+1) - W_O(:,i))
        end do
    end subroutine flux_LF_syst


    subroutine flux_GD_syst (Ns, Flux, W_O)
    ! FLUX POUR SCHEMA DE GODUNOV
        integer, intent(in) :: Ns
        real(rp), dimension(2,Ns), intent(inout) :: Flux
        real(rp), dimension(2,Ns), intent(in) :: W_O
        real(rp) :: v
        integer :: i,j


        do i = 1,(Ns-1)
        !    if (U_O(i) > U_O(i+1)) then ! cas d'une detente pour f concave
        !        if (a_f(U_O(i)) > 0._rp) then
        !            Flux(i) = f(U_O(i))
        !        else if (a_f(U_O(i+1)) < 0._rp) then
        !            Flux(i) = f(U_O(i+1))
        !        else
        !            Flux(i) = a_inv(0._rp)
        !    end if
            !else
            !    v = (f(U_O(i+1)) - f(U_O(i))) / (U_O(i+1) - U_O(i))
            !    if (v > 0._rp) then
            !        Flux(i) = f(U_O(i))
            !    else
            !        Flux(i) = f(U_O(i+1))
            !    end if
            !end if 
        end do

    end subroutine flux_GD_syst

    subroutine flux_RS_syst(Ns, Flux, W_O)
    ! FLUX POUR LE SCHEMA DE RUSANOV
        integer, intent(in) :: Ns
        real(rp), dimension(2,Ns), intent(inout) :: Flux
        real(rp), dimension(2,Ns), intent(in) :: W_O
        real(rp) :: Delta
        integer :: i

        do i = 1,(Ns-1)
            ! max des valeurs propres
            Delta = abs(lambda_1(W_O(:,i))) ! lambda_1 de W_L
            Delta = max(Delta, abs(lambda_1(W_O(:,i+1)))) ! lambda_1 de W_R
            Delta = max(Delta, abs(lambda_2(W_O(:,i)))) ! lambda_2 de W_L
            Delta = max(Delta, abs(lambda_2(W_O(:,i+1)))) ! lambda_2 de W_R
            Flux(:,i) = 0.5_rp*(F(W_O(:,i)) + F(W_O(:,i+1))) - 0.5_rp*Delta * (W_O(:,i+1) - W_O(:,i))
        end do
    end subroutine flux_RS_syst


    subroutine flux_HLL_syst(Ns, Flux, W_O)
        ! FLUX POUR SCHEMA HLL
            integer, intent(in) :: Ns
            real(rp), dimension(2,Ns), intent(inout) :: Flux
            real(rp), dimension(2,Ns), intent(in) :: W_O
            integer :: i
            real(rp) :: bl, br

            do i=1,Ns-1
                ! on calcule b_l et b_r
                bl = min(lambda_1_cons(W_O(:,i)),lambda_1_cons(W_O(:,i+1))) ! min entre lambda_1(W_L) et lambda_1(W_R)
                bl = min(bl,lambda_2_cons(W_O(:,i)),lambda_2_cons(W_O(:,i+1))) ! min entre lambda_2(W_L) et lambda_2(W_R)
                br = max(lambda_1_cons(W_O(:,i)),lambda_1_cons(W_O(:,i+1))) ! max entre lambda_1(W_L) et lambda_1(W_R)
                br = max(br,lambda_2_cons(W_O(:,i)),lambda_2_cons(W_O(:,i+1))) ! max entre lambda_2(W_L) et lambda_2(W_R)
                if (bl > 0.) then
                    Flux(:,i) = F(W_O(:,i))
                else if (bl<0. .AND. br>0.) then
                    Flux(:,i) = (br*W_O(:,i) - bl*W_O(:,i+1) + bl*br*(W_O(:,i+1)-W_O(:,i)))/(br-bl)
                else if (br < 0.) then
                    Flux(:,i) = F(W_O(:,i+1))
                else
                    write(6,*) 'Probleme de calcul de coef pour flux HLL'
                end if
            end do
        end subroutine flux_HLL_syst

        subroutine flux_MR_syst(Ns, Flux, W_O, v_max, rho_max)
            ! FLUX POUR SCHEMA DE MURMAN-ROE
                integer, intent(in) :: Ns
                real(rp), dimension(2,Ns), intent(inout) :: Flux
                real(rp), dimension(2,Ns), intent(in) :: W_O
                real(rp), intent(in) :: v_max, rho_max
                real(rp) :: c
                integer :: i,j


                !do i = 1,(Ns-1)
!
                !    j = 1
                !    if (W_O(j,i) == W_O(j,i+1)) then
                !        c = vitesse(W_O, v_max, rho_max)
                !        c = c*(vitesse(W_O, v_max, rho_max)-W_O(1,i)*p_prime(W_O(1,i),v_max,rho_max))
                !    else
                !        c = (f(W_O(j,i+1))- f(W_O(j,i)))/(W_O(j,i+1) - W_O(j,i))
                !    end if

                !    j = 2
                !    if (W_O(j,i) == W_O(j,i+1)) then
                !        c = vitesse(W_O, v_max, rho_max)
                !        c = c*(vitesse(W_O, v_max, rho_max)-W_O(1,i)*p_prime(W_O(1,i),v_max,rho_max))
                !    else
                !        c = (f(W_O(j,i+1))- f(W_O(j,i)))/(W_O(j,i+1) - W_O(j,i))
                !    end if
                !    Flux(j,i) = 0.5_rp*(f(W_O(i)) + f(W_O(i+1))) - 0.5_rp*abs(c)*(W_O(i+1) - W_O(i))
                !end do
                !end do
            end subroutine flux_MR_syst


end module schemasSW