!
!     surging code in simple & restricted riser surge tank 
!
!     original code : surge.for @ hydraulic foumula progm s3-1
!     ���������W���v���O�����W����13�N�Ły���3-1�z������
!
!*****************************************************************************
    module common_variables
      implicit none
      integer, parameter :: k = kind(1.0d0)
!
      real(k), parameter :: gra = 9.8                   ! �d�͉����x
      real(k), parameter :: pi = 4 * atan(1.)           ! �~���� 3.141592654
!
      real(k), parameter :: q0 = 25.0                   ! �}�ϑO����
      real(k), parameter :: qe = 0.                     ! �}�ό㗬��
      real(k), parameter :: tv = 5.0                    ! �}�ώ���
!
      real(k), parameter :: alp = 1000.0                ! �����H����
      real(k), parameter :: dp = 2.5                    ! �����H���a
      real(k), parameter :: fp = 0.01                   ! �����H���C�����W�� f
      real(k), parameter :: ake = 0.2                   ! �����H���������W�� fe
!
      real(k), parameter :: ds = 7.5                    ! �T�[�W�^���N���a
!
      real(k), parameter :: dth = 1.5                   ! ���������a
      real(k), parameter :: cdp = 0.95                  ! ���������ʌW���i�����F�����H���T�[�W�^���N�j

      real(k), parameter :: tend = 500.                 ! �v�Z����
      integer, parameter :: nt = 1000                   ! ���ԕ�����
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
      ap = pi * dp * dp / 4.                   ! �����H�f�ʐ�
      ak = ake + fp * alp / dp                 ! �����H�����W��
!
      as = pi * ds * ds / 4.                   ! �T�[�W�^���N�f�ʐ�
      ath = pi * dth * dth / 4.                ! �������f�ʐ�
      fath =1.d0/(2.d0*(cdp*ath)**2.d0)        ! �������W��
!
      dt = tend / nt                           ! ���ԍ���
!
      t0 = 2 * pi * sqrt(alp * as / gra / ap)  ! ���R�T�[�W����

      write(1,100) '�����H����L', alp
      write(1,100) '�����H���ad', dp
      write(1,100) '�����H���C�����W��ff', fp
      write(1,100) '�����H���������W��fe', ake
      write(1,100) '�����H�����W��c', ak
      write(1,100) '���R�T�[�W����T0',  t0
      write(1,100)
      write(1,100) '�T�[�W�^���N���aD', ds
      write(1,100) '���������aDp', dth
      write(1,100) '���������ʌW��Cd', cdp
      write(1,100)
      write(1,100) '�}�ϑO����Q0', q0
      write(1,100) '�}�ό㗬��QE', qe
      write(1,100) '�}�ώ���Tv', tv
      write(1,100)
      write(1,100) '�v�Z����Tend', tend
      write(1,110) '���ԕ�����NT', nt
      write(1,100) '���ԍ���DT', dt
      write(1,100)
!
  100 format(a, ':', t25, f12.3)
  110 format(a, ':', t25, i12)
!
      call comp       ! �����Q�E�N�b�^�@�ɂ��v�Z
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
!  �����Q�E�N�b�^�@�ɂ��v�Z
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
      write(1,100) 'T', 'Z', 'V', 'Hp', 'Vp'   ! ����, �T�[�W�^���N����, �����H����, �����������, �������ʉߗ���
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
!  ��b�������̌v�Z
!
!*****************************************************************************
    subroutine func(ff, gg, zz, vv, qq)
      use common_variables
      implicit none
      real(k), intent(in) :: zz, vv, qq
      real(k), intent(out) :: ff, gg
!
      ff = (qq - ap * vv) / as
!      gg = (gra * zz - ak * abs(vv) * v / 2) / alp                                               ! �P���^�T�[�W�^���N
      gg = (gra * zz - ak * abs(vv) * vv / 2 - fath * (ap * vv - qq) * abs(ap * vv - qq)) / alp   ! �������^�T�[�W�^���N
!
      return
!
    end subroutine func
!