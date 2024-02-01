module schemasSH
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
        !call F(F_ex1, W_O(:,1))
        do i = 1,(Ns-1)
            !call F(F_ex, W_O(:,(i+1)))
            !Flux(1,i) = 0.5_rp*(F_ex(1) + F_ex1(1)) - Delta * (W_O(1,i+1) - W_O(1,i))
            !Flux(2,i) = 0.5_rp*(F_ex(2) + F_ex1(2)) - Delta * (W_O(2,i+1) - W_O(2,i))
            Flux(:,i) = 0.5_rp*(F(W_O(:,i)) + F(W_O(:,i+1))) - Delta * (W_O(:,i+1) - W_O(:,i))
            !F_ex1(:) = F_ex(:)
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
        real(rp), dimension(2) :: F_ex1, F_ex
        real(rp) :: Delta
        integer :: i

        !call F(F_ex1, W_O(:,1))
        do i = 1,(Ns-1)
            !call F(F_ex, W_O(:,(i+1)))
            ! max des valeurs propres
            Delta = abs(lambda_1(W_O(:,i))) ! lambda_1 de W_L
            Delta = max(Delta, abs(lambda_1(W_O(:,i+1)))) ! lambda_1 de W_R
            Delta = max(Delta, abs(lambda_2(W_O(:,i)))) ! lambda_2 de W_L
            Delta = max(Delta, abs(lambda_2(W_O(:,i+1)))) ! lambda_2 de W_R
            !Flux(1,i) = 0.5_rp*(F_ex(1) + F_ex1(1)) - 0.5_rp*Delta * (W_O(1,i+1) - W_O(1,i))
            !Flux(2,i) = 0.5_rp*(F_ex(2) + F_ex1(2)) - 0.5_rp*Delta * (W_O(2,i+1) - W_O(2,i))
            !F_ex1(:) = F_ex(:)
            Flux(:,i) = 0.5_rp*(F(W_O(:,i)) + F(W_O(:,i+1))) - 0.5_rp*Delta * (W_O(:,i+1) - W_O(:,i))
        end do
    end subroutine flux_RS_syst


    subroutine flux_HLL_syst(Ns, Flux, W_O)
        ! FLUX POUR SCHEMA HLL
            integer, intent(in) :: Ns
            real(rp), dimension(2,Ns), intent(inout) :: Flux
            real(rp), dimension(2,Ns), intent(in) :: W_O
            real(rp), dimension(2) :: F_ex1, F_ex
            integer :: i
            real(rp) :: l
            !real(rp) :: bl, br

            call F(F_ex1, W_O(:,1))
            do i=1,Ns-1
                ! on calcule b_l et b_r
                !bl = min(lambda_1_cons(W_O(:,i)),lambda_1_cons(W_O(:,i+1)))
                !bl = min(bl,lambda_2_cons(W_O(:,i)),lambda_2_cons(W_O(:,i+1)))
                !br = max(lambda_1_cons(W_O(:,i)),lambda_1_cons(W_O(:,i+1)))
                !br = max(br,lambda_2_cons(W_O(:,i)),lambda_2_cons(W_O(:,i+1)))
                !if (bl >= 0) then
                !    call F(F_ex, W_O(:,(i)))
                !    Flux(:,i) = F_ex(:)
                !else if (bl<0 .AND. br>=0) then
                !    call F(F_ex, W_O(:,(i)))
                !    call F(F_ex1, W_O(:,(i+1)))
                !    Flux(1,i) = (br*F_ex1(1)-bl*F_ex(1)+(bl*br)*(W_O(1,i+1)-W_O(1,i)))/(br-bl)
                !    Flux(2,i) = (br*F_ex1(2)-bl*F_ex(2)+(bl*br)*(W_O(2,i+1)-W_O(2,i)))/(br-bl)
                !else if (br < 0) then
                !    call F(F_ex, W_O(:,(i+1)))
                !    Flux(:,i) = F_ex(:)
                !else
                !    write(6,*) 'Probleme de calcul de coef pour flux HLL'
                !end if
                l = max(lambda_2_cons(W_O(:,i)),lambda_2_cons(W_O(:,i+1)))  
                Flux(:,i) = 0.5_rp*(F(W_O(:,i)) + F(W_O(:,i+1))) - 0.5_rp*l*(W_O(:,i)-W_O(:,i+1))
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


end module schemasSH