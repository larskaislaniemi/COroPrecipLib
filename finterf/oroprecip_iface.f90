module m_oroprecip
    implicit none

    type :: oroprecip_options
        integer :: nx, ny
        double precision :: dx, dy, Lx, Ly
        double precision :: U(2)
        double precision :: tauc, tauf, Sbg
        integer :: sink_at_downslope
    end type oroprecip_options

    type oroprecip_thermodyn
        double precision :: rho_ref, qsat_ref, Hm
    end type oroprecip_thermodyn

    interface
        subroutine oroprecip(h, co, ct, p)
            use, intrinsic :: iso_c_binding
            import :: oroprecip_options
            import :: oroprecip_thermodyn

            implicit none

            real (kind=c_double), dimension(*), intent(in) :: h
            real (kind=c_double), dimension(*), intent(out) :: p
            type (oroprecip_options), intent(in) :: co
            type (oroprecip_thermodyn), intent(in) :: ct
        end subroutine oroprecip
    end interface
end module m_oroprecip

