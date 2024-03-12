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
        real(rp) :: Delta
        integer :: i

        ! calcul flux pour Lax-Friedrichs
        Delta = (dx / dt) * 0.5_rp
        do i = 1,(Ns-1)
            Flux(:,i) = 0.5_rp*(F(W_O(:,i)) + F(W_O(:,i+1))) - Delta * (W_O(:,i+1) - W_O(:,i))
        end do
    end subroutine flux_LF_syst


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


    subroutine flux_HLL_syst(Ns, Flux, W_O, dx, dt, lambda)
        ! FLUX POUR SCHEMA HLL
        integer, intent(in) :: Ns
        real(rp), intent(in) :: dt, dx
        real(rp), dimension(2,Ns), intent(inout) :: Flux
        real(rp), dimension(2,Ns), intent(in) :: W_O
        real(rp), intent(out) :: lambda
        integer :: i
        real(rp) :: hl, hr, ul, ur, pil, pir

        lambda = 2.*dt/dx
        do i=1,Ns-1
            hl = W_O(1,i)
            hr = W_O(1,i+1)
       
            ul = vitesse(W_O(:,i))
            ur = vitesse(W_O(:,i+1))
       
            pil = hl*ul**2 + g*0.5*hl**2
            pir = hr*ur**2 + g*0.5*hr**2
       
            Flux(1,i) = 0.5*(hl*ul+hr*ur) - 0.5/lambda*(hr-hl)
            Flux(2,i) = 0.5*(pil + pir) -0.5/lambda*(hr*ur-hl*ul)
        end do
    end subroutine flux_HLL_syst


    subroutine flux_recons_hydro(Ns, Flux, W_O, Zi, dx, dt, W_Om, W_Op)
        ! FLUX POUR LA RECONSTRUCTION HYDROSTATIQUE
        integer, intent(in) :: Ns
        real(rp), intent(in) :: dt, dx
        real(rp), dimension(2,Ns), intent(inout) :: Flux
        real(rp), dimension(2,Ns), intent(in) :: W_O
        real(rp), dimension(Ns), intent(in) :: Zi
        real(rp), dimension(2,Ns), intent(out) :: W_Op
        real(rp), dimension(2,Ns), intent(out) :: W_Om
        integer :: i
        real(rp) :: z_idemi, h_idemip, h_idemim, lambda
        real(rp) :: hl, hr, ul, ur, pil, pir

        do i = 1,Ns-1
            z_idemi = max(Zi(i), Zi(i+1)) ! on definit z_{i+1/2}
            ! on definit h_{i+1/2}^+ et h_{i+1/2}^-
            h_idemim = max(0._rp, W_O(1,i)+(Zi(i)-z_idemi))
            h_idemip = max(0._rp, W_O(1,i+1)+(Zi(i+1)-z_idemi))
            W_Om(1,i) = h_idemim
            W_Om(2,i) = h_idemim*(W_O(2,i)/W_O(1,i))
            W_Op(1,i) = h_idemip
            W_Op(2,i) = h_idemip*(W_O(2,i+1)/W_O(1,i+1))
        end do

        lambda = 2.*dt/dx
        ! flux HLL
        do i = 1,Ns-1
            hl = W_Om(1,i)
            hr = W_Op(1,i)

            ul = vitesse(W_Om(:,i))
            ur = vitesse(W_Op(:,i))

            pil = hl*ul**2 + g*0.5*hl**2
            pir = hr*ur**2 + g*0.5*hr**2

            Flux(1,i) = 0.5*(hl*ul+hr*ur) - 0.5/lambda*(hr-hl)
            Flux(2,i) = 0.5*(pil + pir) -0.5/lambda*(hr*ur-hl*ul)
        end do

    end subroutine flux_recons_hydro

    subroutine flux_GDWB(Ns, Flux, W_O, Zi, dx, dt)
        ! FLUX POUR LE SCHEMA TYPE GODUNOV WB
        integer, intent(in) :: Ns
        real(rp), intent(in) :: dt, dx
        real(rp), dimension(2,Ns), intent(inout) :: Flux
        real(rp), dimension(2,Ns), intent(in) :: W_O
        real(rp), dimension(Ns), intent(in) :: Zi
        integer :: i
        real(rp) :: lambda, pil,pir,hl,hr,ul,ur,zr,zl

        lambda = 2.*dt/dx
        do i =1,Ns-1
            hl = W_O(1,i)
            hr = W_O(1,i+1)
       
            ul = vitesse(W_O(:,i))
            ur = vitesse(W_O(:,i+1))
       
            pil = hl*ul**2 + g*0.5*hl**2
            pir = hr*ur**2 + g*0.5*hr**2
            zr = Zi(i+1)
            zl = Zi(i)
       
            Flux(1,i) = 0.5*(hl*ul+hr*ur) - 0.5/lambda*((hr+zr)-(hl+zl))
            Flux(2,i) = 0.5*(pil + pir) -0.5/lambda*(hr*ur-hl*ul)
        end do
    end subroutine flux_GDWB

end module schemasSW