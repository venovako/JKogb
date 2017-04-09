module h5data
  use h5lt
  use hdf5
  implicit none

  integer, parameter :: idalen = 4
  integer, parameter :: timlen = 4

contains

  subroutine readh0(gid, idadim, info)
    implicit none

    integer(hid_t), intent(in) :: gid
    integer, intent(out) :: idadim(idalen), info(2)

    integer(hsize_t) :: dims1(1)

    info = 0

    dims1(1) = idalen
1   call h5ltread_dataset_int_f(gid, 'IDADIM', idadim, dims1, info(2))
    if (info(2) .ne. 0) info(1) = 1
  end subroutine readh0

  subroutine readh1(gid, m, n, G, ldg, jvec, info)
    implicit none

    integer(hid_t), intent(in) :: gid
    integer, intent(in) :: m, n, ldg
    double precision, intent(out) :: G(ldg, n)
    integer, intent(out) :: jvec(n), info(2)

    integer(hsize_t) :: dims1(1), dims2(2)

    info = 0

    dims2(1) = ldg
    dims2(2) = n
1   call h5ltread_dataset_double_f(gid, 'G', G, dims2, info(2))
    if (info(2) .ne. 0) then
       info(1) = 1
       return
    end if

    dims1(1) = n
2   call h5ltread_dataset_int_f(gid, 'JVEC', jvec, dims1, info(2))
    if (info(2) .ne. 0) then
       info(1) = 2
       return
    end if
  end subroutine readh1

  subroutine readh2(gid, n, ipl, invp, info)
    implicit none

    integer(hid_t), intent(in) :: gid
    integer, intent(in) :: n
    integer, intent(out) :: ipl(n), invp(n), info(2)

    integer(hsize_t) :: dims1(1)

    info = 0

    dims1(1) = n
1   call h5ltread_dataset_int_f(gid, 'IPL', ipl, dims1, info(2))
    if (info(2) .ne. 0) then
       info(1) = 1
       return
    end if

    dims1(1) = n
2   call h5ltread_dataset_int_f(gid, 'INVP', invp, dims1, info(2))
    if (info(2) .ne. 0) then
       info(1) = 2
       return
    end if
  end subroutine readh2

  subroutine dumph5(gid, m, n, G, ldg, jvec, nplus, U, ldu, V, ldv, Q, ldq, tol, iflags, status, Tcpu, info)
    implicit none

    integer(hid_t), intent(in) :: gid
    integer, intent(in) :: m, n, ldg, jvec(n), nplus, ldu, ldv, ldq, iflags, status
    double precision, intent(in) :: G(ldg, n), U(ldu, m), V(ldv, n), Q(ldq, m), tol, Tcpu(timlen)
    integer, intent(out) :: info(2)

    integer(hsize_t) :: dims1(1), dims2(2)
    integer(size_t) :: atrlen

    integer :: idummy(3)
    double precision :: ddummy(1)

    info = 0

    if (ldg .gt. 0) then
       dims2(1) = ldg
       dims2(2) = n
1      call h5ltmake_dataset_double_f(gid, 'G', 2, dims2, G, info(2))
       if (info(2) .ne. 0) then
          info(1) =  1
          return
       end if

       atrlen = 3
       idummy(1) = m
       idummy(2) = n
       idummy(3) = ldg
2      call h5ltset_attribute_int_f(gid, 'G', 'dims', idummy, atrlen, info(2))
       if (info(2) .ne. 0) then
          info(1) =  2
          return
       end if

       dims1(1) = n
3      call h5ltmake_dataset_int_f(gid, 'JVEC', 1, dims1, jvec, info(2))
       if (info(2) .ne. 0) then
          info(1) =  3
          return
       end if

       atrlen = 3
       idummy(1) = n
       idummy(2) = nplus
       idummy(3) = m
4      call h5ltset_attribute_int_f(gid, 'JVEC', 'dims', idummy, atrlen, info(2))
       if (info(2) .ne. 0) then
          info(1) =  4
          return
       end if
    end if

    if (ldu .gt. 0) then
       dims2(1) = ldu
       dims2(2) = m
5      call h5ltmake_dataset_double_f(gid, 'U', 2, dims2, U, info(2))
       if (info(2) .ne. 0) then
          info(1) =  5
          return
       end if

       atrlen = 3
       idummy(1) = m
       idummy(2) = m
       idummy(3) = ldu
6      call h5ltset_attribute_int_f(gid, 'U', 'dims', idummy, atrlen, info(2))
       if (info(2) .ne. 0) then
          info(1) =  6
          return
       end if
    end if

    if (ldv .gt. 0) then
       dims2(1) = ldv
       dims2(2) = n
7      call h5ltmake_dataset_double_f(gid, 'V', 2, dims2, V, info(2))
       if (info(2) .ne. 0) then
          info(1) =  7
          return
       end if

       atrlen = 3
       idummy(1) = n
       idummy(2) = n
       idummy(3) = ldv
8      call h5ltset_attribute_int_f(gid, 'V', 'dims', idummy, atrlen, info(2))
       if (info(2) .ne. 0) then
          info(1) =  8
          return
       end if
    end if

    if (ldq .gt. 0) then
       dims2(1) = ldq
       dims2(2) = m
9      call h5ltmake_dataset_double_f(gid, 'Q', 2, dims2, Q, info(2))
       if (info(2) .ne. 0) then
          info(1) =  9
          return
       end if

       atrlen = 3
       idummy(1) = m
       idummy(2) = m
       idummy(3) = ldq
10     call h5ltset_attribute_int_f(gid, 'Q', 'dims', idummy, atrlen, info(2))
       if (info(2) .ne. 0) then
          info(1) = 10
          return
       end if
    end if

    dims1(1) = timlen
11  call h5ltmake_dataset_double_f(gid, 'Tcpu', 1, dims1, Tcpu, info(2))
    if (info(2) .ne. 0) then
       info(1) = 11
       return
    end if

    atrlen = 1
    idummy(1) = timlen
12  call h5ltset_attribute_int_f(gid, 'Tcpu', 'dims', idummy, atrlen, info(2))
    if (info(2) .ne. 0) then
       info(1) = 12
       return
    end if

    dims1(1) = 2
    idummy(1) = status
    idummy(2) = iflags
13  call h5ltmake_dataset_int_f(gid, 'STATUS', 1, dims1, idummy, info(2))
    if (info(2) .ne. 0) then
       info(1) = 13
       return
    end if

    atrlen = 1
    ddummy(1) = tol
14  call h5ltset_attribute_double_f(gid, 'STATUS', 'tol', ddummy, atrlen, info(2))
    if (info(2) .ne. 0) then
       info(1) = 14
       return
    end if

    atrlen = 1
    idummy(1) = 2
15  call h5ltset_attribute_int_f(gid, 'STATUS', 'dims', idummy, atrlen, info(2))
    if (info(2) .ne. 0) then
       info(1) = 15
       return
    end if
  end subroutine dumph5
end module h5data
