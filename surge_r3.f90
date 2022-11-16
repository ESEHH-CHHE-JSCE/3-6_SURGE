!
!     surging code in simple & restricted riser surge tank 
!
!     original code : surge.for @ hydraulic foumula progm s3-1
!     …—Œö®W—á‘èƒvƒƒOƒ‰ƒ€W•½¬13”N”Åy—á‘è3-1z‚ğ‰ü•Ï
!
!*****************************************************************************
    module common_variables
      implicit none
      integer, parameter :: k = kind(1.0d0)
!
      real(k), parameter :: gra = 9.8                   ! d—Í‰Á‘¬“x
      real(k), parameter :: pi = 4 * atan(1.)           ! ‰~ü—¦ 3.141592654
!
      real(k), parameter :: q0 = 25.0                   ! ‹}•Ï‘O—¬—Ê
      real(k), parameter :: qe = 0.                     ! ‹}•ÏŒã—¬—Ê
      real(k), parameter :: tv = 5.0                    ! ‹}•ÏŠÔ
!
      real(k), parameter :: alp = 1000.0                ! “±…˜H’·‚³
      real(k), parameter :: dp = 2.5                    ! “±…˜H“àŒa
      real(k), parameter :: fp = 0.01                   ! “±…˜H–€C‘¹¸ŒW” f
      real(k), parameter :: ake = 0.2                   ! “±…˜H“üŒû‘¹¸ŒW” fe
!
      real(k), parameter :: ds = 7.5                    ! ƒT[ƒWƒ^ƒ“ƒN’¼Œa
!
      real(k), parameter :: dth = 1.5                   ! §…Œû’¼Œa
      real(k), parameter :: cdp = 0.95                  ! §…Œû—¬—ÊŒW”i—¬“üF“±…˜H¨ƒT[ƒWƒ^ƒ“ƒNj

      real(k), parameter :: tend = 500.                 ! ŒvZŠÔ
      integer, parameter :: nt = 1000                   ! ŠÔ•ªŠ„”
!
      real(k) :: ap, ak, as, ath, fath, dt
!
    end module common_variables
!
!*****************************************************************************
    program surge
      use common_variables
      implicit none
      real(k) :: t0
!
      open(1,file='surge.dat',status='unknown')
!
      ap = pi * dp * dp / 4.                   ! “±…˜H’f–ÊÏ
      ak = ake + fp * alp / dp                 ! “±…˜H‘¹¸ŒW”
!
      as = pi * ds * ds / 4.                   ! ƒT[ƒWƒ^ƒ“ƒN’f–ÊÏ
      ath = pi * dth * dth / 4.                ! §…Œû’f–ÊÏ
      fath =1.d0/(2.d0*(cdp*ath)**2.d0)        ! §…ŒûŒW”
!
      dt = tend / nt                           ! ŠÔ‚İ
!
      t0 = 2 * pi * sqrt(alp * as / gra / ap)  ! ©—RƒT[ƒWüŠú

      write(1,100) '“±…˜H’·‚³L', alp
      write(1,100) '“±…˜H“àŒad', dp
      write(1,100) '“±…˜H–€C‘¹¸ŒW”ff', fp
      write(1,100) '“±…˜H—¬“ü‘¹¸ŒW”fe', ake
      write(1,100) '“±…˜H‘¹¸ŒW”c', ak
      write(1,100) '©—RƒT[ƒWüŠúT0',  t0
      write(1,100)
      write(1,100) 'ƒT[ƒWƒ^ƒ“ƒN’¼ŒaD', ds
      write(1,100) '§…Œû’¼ŒaDp', dth
      write(1,100) '§…Œû—¬—ÊŒW”Cd', cdp
      write(1,100)
      write(1,100) '‹}•Ï‘O—¬—ÊQ0', q0
      write(1,100) '‹}•ÏŒã—¬—ÊQE', qe
      write(1,100) '‹}•ÏŠÔTv', tv
      write(1,100)
      write(1,100) 'ŒvZŠÔTend', tend
      write(1,110) 'ŠÔ•ªŠ„”NT', nt
      write(1,100) 'ŠÔ‚İDT', dt
      write(1,100)
!
  100 format(a, ':', t25, f12.3)
  110 format(a, ':', t25, i12)
!
      call comp       ! ƒ‹ƒ“ƒQEƒNƒbƒ^–@‚É‚æ‚éŒvZ
!
      close(1)
!
      stop
!
    end program surge
!
!*****************************************************************************
!
!  subroutine: comp
!
!  ƒ‹ƒ“ƒQEƒNƒbƒ^–@‚É‚æ‚éŒvZ
!
!*****************************************************************************
    subroutine comp
      use common_variables
      implicit none
      real(k) :: f1, f2, f3, f4
      real(k) :: g1, g2, g3, g4
      real(k) :: z1, z2, z3, z4
      real(k) :: v1, v2, v3, v4
      real(k) :: qq, tt, v0, z0, zn, vn
      integer :: i
!
      write(1,100) 'T', 'Z', 'V', 'Hp', 'Vp'   ! ŠÔ, ƒT[ƒWƒ^ƒ“ƒN…ˆÊ, “±…˜H—¬‘¬, §…ŒûŠî•”…“ª, §…Œû’Ê‰ß—¬‘¬
!
      tt = 0.
      v0 = q0 / ap
!
      z0 = ak * v0 ** 2 / 2 / gra
      write(1,110) tt, z0, v0, z0, 0
!
      do 10 i = 0, nt - 1
          if(tt.le.tv) then
              qq = q0 + (qe - q0)/tv*tt
          else
              qq = qe
          endif
!
          call func(f1, g1, z0, v0, qq)
          z1 = z0 + f1 * dt / 2
          v1 = v0 + g1 * dt / 2
!
          call func(f2, g2, z1, v1, qq)
          z2 = z0 + f2 * dt / 2
          v2 = v0 + g2 * dt / 2

          call func(f3, g3, z2, v2, qq)
!          z3 = z0 + f3 * dt / 2     ! surge.for
!          v3 = v0 + g3 * dt / 2     ! surge.for
          z3 = z0 + f3 * dt
          v3 = v0 + g3 * dt
!
          call func(f4, g4, z3, v3, qq)
!
          zn = z0 + (f1 + 2 * f2 + 2 * f3 + f4) * dt / 6
          vn = v0 + (g1 + 2 * g2 + 2 * g3 + g4) * dt / 6
!
          tt = (i + 1) * dt
          write(1,110) tt, zn, vn, zn - ((ap * vn - qq) * abs(ap * vn - qq)/(ath**2.)/2./gra), (ap * vn - qq)
!
          z0 = zn
          v0 = vn
!
   10 continue
!
  100 format(t8,a,t20,a,t32,a,t44,a,t56,a)
  110 format(5f12.3)
!
      return
!
    end subroutine comp
!
!*****************************************************************************
!
!  subroutine: func
!
!  Šî‘b•û’ö®‚ÌŒvZ
!
!*****************************************************************************
    subroutine func(ff, gg, zz, vv, qq)
      use common_variables
      implicit none
      real(k), intent(in) :: zz, vv, qq
      real(k), intent(out) :: ff, gg
!
      ff = (qq - ap * vv) / as
!      gg = (gra * zz - ak * abs(vv) * v / 2) / alp                                               ! ’P“®Œ^ƒT[ƒWƒ^ƒ“ƒN
      gg = (gra * zz - ak * abs(vv) * vv / 2 - fath * (ap * vv - qq) * abs(ap * vv - qq)) / alp   ! §…ŒûŒ^ƒT[ƒWƒ^ƒ“ƒN
!
      return
!
    end subroutine func
!