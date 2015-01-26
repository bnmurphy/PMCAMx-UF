MODULE io_ezcdf
  !!
  !! Netcdf input/output
  !! Still in fortran77 netcdf include style, but works...
  !! 
  !!
  !!
  !! [ eZcdf ]
  !!
  !! Last update, Brodeau, 2009
  !!
  !!
  !!
  IMPLICIT NONE
  !!
  PRIVATE
  !!
  !! List of available routines
  !! ==========================
  PUBLIC :: dims,             &
       &    get_sf_ao,        &
       &    getvar_1d,        &
       &    getvar_2d,        &
       &    getvar_3d,        &
       &    getmask,          &
       &    getmask_3d,       &
       &    p2d_t, p2d_t_irr, &
       &    p3d_t, p3d_t_irr, &
       &    check_4_miss,     &
       &    get_var_info,     &
       &    prtmask
  !!
  !!===========================
  !!
  !!
  INCLUDE 'netcdf.inc'
  !!
  CHARACTER(len=20), PUBLIC :: cv_misc
  !!
  INTEGER, PARAMETER, PUBLIC :: &
       &    wpez = 4
  !!
  INTEGER, PUBLIC :: &
       &    ils,           &
       &    n_dims,        &
       &    nvid,        &
       &    nfid,          &   !: netcdf ID file for only 1 snapshot use !!!
       &    id_x, id_y,    &
       &    id_z, id_t,    &
       &    id_lon, id_lat,&
       &    id_dpth, id_time
  !!
  REAL(wpez), DIMENSION(2) :: &
       &    rrange
  !!
  INTEGER, DIMENSION(2) :: &
       &    id_dims2,    &
       &    istart2,     &
       &    icount2
  !!
  INTEGER, DIMENSION(3) :: &
       &    id_dims3,   &
       &    istart3,    &
       &    icount3
  !!
  INTEGER, DIMENSION(4) :: & 
       &    id_dims4,   &
       &    istart4,    &
       &    icount4
  !!
  REAL(wpez)            :: & 
       &    rmin
  !!
  REAL(wpez), PARAMETER :: &
       &    vflagout = -9999.
  !!
  CHARACTER(len=100)     :: &
       &    croutnm,       &
       &    cu
  !!
  CHARACTER(LEN=400), PARAMETER   ::     &
       &    cabout = 'File created by SOSIE interpolation environement,&
       & Laurent Brodeau, 2008, (brodeau@gmail.com)'
  !! 
  CHARACTER(len=20) :: cvxout, cvyout
  !!
  INTEGER :: &
       &  ji, jj
  !!
  !!
CONTAINS
  !!
  !!========================================================================
  !!
  !!
  SUBROUTINE DISP_ERR(ierr, croutine, ctf, ctv)
    !!
    !! To handle and display error messages
    !!
    INTEGER,            INTENT(in) :: ierr
    !!
    CHARACTER(len=*) , INTENT(in) :: &
         &            ctf,           &    !: treated file
         &            croutine,      &    !: routine name
         &            ctv                 !: treated varible
    !!
    !!
    IF ( ierr /= 0 ) THEN
       PRINT *, ''
       PRINT *, '************************************************'
       PRINT *, 'Error occured in procedure ', trim(croutine),' !'
       PRINT *, ''
       PRINT *, 'Treated file     = ', trim(ctf)
       PRINT *, 'Treated variable = ', trim(ctv)
       PRINT *, ''
       PRINT *, '--> aborting program'
       PRINT *, ''
       PRINT *, 'Netcdf message was :'
       PRINT *, trim(NF_STRERROR(ierr))
       PRINT *, ''
       PRINT *, '************************************************'
       PRINT *, ''
       STOP
    END IF
    !!
  END SUBROUTINE DISP_ERR
  !!
  !!===============================================================================
  !!
  !!
  SUBROUTINE DIMS(cfil, cvar, lx, ly, lz, lt)
    !!
    !!-----------------------------------------------------------------------
    !! This routine opens a netcdf file 'cfil' to check the dimension 
    !! of the variable 'cvar'. It then gives the length of each of the dimension,
    !! if the length returns '-1' that means that the dimension does not exist
    !!
    !! example : if the variable has only 1 dimension, of length 132,
    !!           DIMS will return lx=132, ly=-1, lz=-1, lt=-1
    !!
    !!
    !! INPUT :
    !! -------
    !!          * cfil       : name of the input file          (character)
    !!          * cvar       : name of the variable            (character)
    !!
    !! OUTPUT :
    !! --------
    !!          * lx      : record length of first dimension       (integer)
    !!          * ly      : record length of second dimension      (integer)
    !!          * lz      : record length of third dimension       (integer)
    !!          * lt      : record length of fourth dimension      (integer)
    !!
    !!
    !!  Author :            Laurent Brodeau, 2004
    !!  --------
    !!------------------------------------------------------------------------
    !!
    CHARACTER(len=*), INTENT(in)  :: cfil, cvar
    INTEGER,            INTENT(out) :: lx, ly, lz, lt           
    !!
    INTEGER :: jdim
    !!
    !!
    croutnm = 'DIMS'
    !!
    CALL disp_err(NF_OPEN(cfil, NF_NOWRITE, nfid), croutnm, cfil, cvar)
    !!
    !! Querrying dimensions on variable :
    CALL disp_err(NF_INQ_VARID(nfid,cvar,nvid), croutnm, cfil, cvar)
    CALL disp_err(NF_INQ_VARNDIMS(nfid,nvid,n_dims), croutnm, cfil, cvar)
    !! 
    !!        1 dimension
    IF ( n_dims == 1 ) THEN
       CALL disp_err(NF_INQ_VARDIMID(nfid, nvid, jdim), croutnm, cfil, cvar)
       CALL disp_err(NF_INQ_DIMLEN(nfid, jdim, lx), croutnm, cfil, cvar)
       ly = -1 ; lz = -1 ; lt = -1
    END IF
    !!
    !!        2 dimensions
    IF ( n_dims == 2 ) THEN 
       CALL disp_err(NF_INQ_VARDIMID(nfid,nvid,id_dims2), croutnm, cfil, cvar)
       !!
       CALL disp_err(NF_INQ_DIMLEN(nfid,id_dims2(1),lx), croutnm, cfil, cvar)  
       CALL disp_err(NF_INQ_DIMLEN(nfid,id_dims2(2),ly), croutnm, cfil, cvar)
       !!
       lz = -1 ; lt = -1
       !!
    END IF
    !!
    !!        3 dimensions
    IF ( n_dims == 3 ) THEN 
       CALL disp_err(NF_INQ_VARDIMID(nfid,nvid,id_dims3), croutnm, cfil, cvar)
       CALL disp_err(NF_INQ_DIMLEN(nfid,id_dims3(1),lx), croutnm, cfil, cvar)  
       CALL disp_err(NF_INQ_DIMLEN(nfid,id_dims3(2),ly), croutnm, cfil, cvar)
       CALL disp_err(NF_INQ_DIMLEN(nfid,id_dims3(3),lz), croutnm, cfil, cvar)     
       lt = -1
    END IF
    !!
    !!        4 dimensions
    IF ( n_dims == 4 ) THEN 
       CALL disp_err(NF_INQ_VARDIMID(nfid,nvid,id_dims4), croutnm, cfil, cvar)
       CALL disp_err(NF_INQ_DIMLEN(nfid,id_dims4(1),lx), croutnm, cfil, cvar)  
       CALL disp_err(NF_INQ_DIMLEN(nfid,id_dims4(2),ly), croutnm, cfil, cvar)
       CALL disp_err(NF_INQ_DIMLEN(nfid,id_dims4(3),lz), croutnm, cfil, cvar)
       CALL disp_err(NF_INQ_DIMLEN(nfid,id_dims4(4),lt), croutnm, cfil, cvar)     
    END IF
    !!
    IF ( n_dims > 4 ) THEN
       PRINT *, 'Problem into eZcdf routine DIMS:'
       PRINT *, ' - the dimension is > 4!', n_dims
       STOP
    END IF
    !!
    IF ( n_dims < 1 ) THEN
       PRINT *, 'Problem into eZcdf routine DIMS: '
       PRINT *, 'The dimension is unexpected!', n_dims
       STOP
    END IF
    !!
    !! Closing file :
    CALL disp_err(NF_CLOSE(nfid), croutnm, cfil, cvar)
    !!
  END SUBROUTINE DIMS
  !!
  !!===============================================================================
  !!
  !!
  SUBROUTINE GETVAR_1D(cfil, cvar, lx, X)
    !!
    !!-----------------------------------------------------------------------
    !! This routine extract a variable 1D from a netcdf file
    !!
    !! INPUT :
    !! -------
    !!          * cfil      : name of the input file             (character l=100)
    !!          * cvar      : name of the variable               (character l=20)
    !!
    !! OUTPUT :
    !! --------
    !!          * lx        : length dimension of the variable  (integer)
    !!          * X         : 1D array contening the variable    (real)
    !!
    !!
    !!  Author :            Laurent Brodeau, 2004
    !!  --------
    !!------------------------------------------------------------------------
    !!
    CHARACTER(len=*), INTENT(in) :: cfil, cvar
    INTEGER,            INTENT(in) :: lx
    REAL(wpez), DIMENSION (lx), INTENT(out) ::  X
    !!
    croutnm = 'GETVAR_1D'
    !!
    CALL disp_err(NF_OPEN(cfil, NF_NOWRITE, nfid), croutnm, cfil, cvar)
    CALL disp_err(NF_INQ_VARID(nfid, trim(cvar), nvid), croutnm, cfil, cvar)
    CALL disp_err(NF_GET_VAR_REAL(nfid,nvid,X), croutnm, cfil, cvar)
    !!        
    !! Closing file :
    CALL disp_err(NF_CLOSE(nfid), croutnm, cfil, cvar)
    !!
  END SUBROUTINE GETVAR_1D
  !!
  !!===============================================================================
  !!
  !!
  SUBROUTINE GETVAR_2D(id_fil, id_var, cfil, cvar, lx, ly, lt, kz, kt, X, jt1, jt2)
    !!
    !!-----------------------------------------------------------------------------
    !! This routine extract a variable 2D from a netcdf file
    !! at a given time
    !!
    !! INPUT :
    !! -------
    !!          * id_fil    : ID of current file                  (integer)
    !!          * id_var    : ID of current variable              (integer)
    !!          * cfil      : name of the input file              (character)
    !!          * cvar      : name of the variable                (character)
    !!          * lx        : x dimension of the variable        (integer)
    !!          * ly        : y dimension of the variable        (integer)
    !!          * lt        : time dimension of the variable     (integer)
    !!          * kz        : level to extract                    (integer)
    !!                      0 => input file does not have levels
    !!
    !!          * kt        : time snapshot to extract            (integer)
    !!                      0 => input file does not have a time snapshot
    !!                           (= old GETVAR_2D_NOTIME)
    !!
    !!
    !! OUTPUT :
    !! --------
    !!          * X         : 2D array contening the variable     (real)
    !!
    !! OPTIONAL INPUT :
    !! ----------------
    !!          * jt1, jt2  : first and last time snapshot to extract
    !!
    !!
    !!  Author :            Laurent Brodeau, 2004
    !!  --------
    !!------------------------------------------------------------------------
    !!
    INTEGER,                     INTENT(inout)  :: id_fil, id_var
    CHARACTER(len=*),             INTENT(in)  :: cfil, cvar
    INTEGER,                        INTENT(in)  :: lx, ly, lt, kz, kt
    REAL(wpez), DIMENSION (lx, ly), INTENT(out) :: X
    INTEGER,              OPTIONAL, INTENT(in)  :: jt1, jt2
    !!
    INTEGER :: its, ite
    !!
    croutnm = 'GETVAR_2D'
    !!
    !!
    IF ( present(jt1).AND.present(jt2) ) THEN
       its = jt1 ; ite = jt2
    ELSE
       its = 1   ; ite = lt
    END IF
    !!
    !!
    IF ( (kt == its).OR.(kt == 0) ) THEN   ! Opening file and defining variable :
       CALL disp_err(NF_OPEN(cfil, NF_NOWRITE, id_fil), croutnm, cfil, cvar)
       CALL disp_err(NF_INQ_VARID(id_fil, cvar, id_var), croutnm, cfil, cvar)
    END IF
    !!
    !!
    IF ( kz == 0 ) THEN    ! No levels
       !!
       !!
       IF ( kt == 0 ) THEN
          !!
          istart2 = (/ 1  , 1  /)
          icount2 = (/ lx , ly /)
          CALL disp_err(NF_GET_VARA_REAL(id_fil,id_var,istart2,icount2,X), &
               &        croutnm, cfil, cvar)
          !!
       ELSEIF ( kt > 0 ) THEN
          !!
          istart3 = (/ 1 , 1 , kt /)
          icount3 = (/ lx , ly , 1 /)
          CALL disp_err(NF_GET_VARA_REAL(id_fil,id_var,istart3,icount3,X), &
               &        croutnm, cfil, cvar)
          !!
       END IF
       !!
       !!
    ELSEIF ( kz > 0 ) THEN
       !!
       !!
       IF ( kt == 0 ) THEN
          !!
          istart3 = (/ 1  , 1  , kz /)
          icount3 = (/ lx , ly , 1  /)
          CALL disp_err(NF_GET_VARA_REAL(id_fil,id_var,istart3,icount3,X), &
               &        croutnm, cfil, cvar)
          !!
       ELSEIF ( kt > 0 ) THEN
          !!
          istart4 = (/ 1 , 1 , kz , kt /)
          icount4 = (/ lx , ly , 1 , 1 /)
          CALL disp_err(NF_GET_VARA_REAL(id_fil,id_var,istart4,icount4,X), &
               &        croutnm, cfil, cvar)          
          !!
       END IF
       !!
       !!
    END IF
    !!
    IF (( kt == ite ).or.( kt == 0 ))  THEN
       CALL disp_err(NF_CLOSE(id_fil), croutnm, cfil, cvar)
       id_fil = 0 
       id_var = 0
    END IF
    !!
  END SUBROUTINE GETVAR_2D
  !!
  !!===============================================================================
  !!
  !!
  !!
  !!
  SUBROUTINE GETVAR_3D(id_fil, id_var, cfil, cvar, lx, ly, lz, lt, kt, X)
    !!
    !!------------------------------------------------------------------
    !! This routine extract a variable 3D from a netcdf file
    !! at a given time
    !!
    !! INPUT :
    !! -------
    !!          * id_fil    : ID of current file                  (integer)
    !!          * id_var    : ID of current variable              (integer)
    !!          * cfil      : name of the input file              (character)
    !!          * cvar      : name of the variable                (character)
    !!          * lx        : x dimension of the variable        (integer)
    !!          * ly        : y dimension of the variable        (integer)
    !!          * lz        : z dimension of the variable        (integer)
    !!          * lt        : time dimension of the variable     (integer)
    !!
    !!          * kt        : time snapshot to extract            (integer)
    !!                      0 => input file does not have a time snapshot
    !!                           (= old GETVAR_2D_NOTIME)
    !!
    !!
    !! OUTPUT :
    !! --------
    !!          * X         : 3D array contening the variable     (real)
    !!
    !!
    !!  Author :            Laurent Brodeau, 2004
    !!  --------
    !!------------------------------------------------------------------------
    !!
    INTEGER,                         INTENT(inout)  :: id_fil, id_var
    CHARACTER(len=*),                 INTENT(in)  :: cfil, cvar
    INTEGER,                            INTENT(in)  :: lx, ly, lz, lt, kt
    REAL(wpez), DIMENSION (lx, ly, lz), INTENT(out) :: X
    !!
    croutnm = 'GETVAR_3D'
    !!
    !!
    IF ( kt <= 1 ) THEN   ! Opening file and defining variable :
       CALL disp_err(NF_OPEN(cfil, NF_NOWRITE,  id_fil), croutnm, cfil, cvar)
       CALL disp_err(NF_INQ_VARID(id_fil, cvar, id_var), croutnm, cfil, cvar)
    END IF
    !!
    !!
    IF ( kt == 0 ) THEN
       !!
       istart3 = (/ 1  , 1  , 1  /)
       icount3 = (/ lx , ly , lz /)
       CALL disp_err(NF_GET_VARA_REAL(id_fil, id_var, istart3, icount3, X), &
            &        croutnm, cfil, cvar)
       !!
    ELSEIF ( kt > 0 ) THEN
       !!
       istart4 = (/ 1 , 1 , 1 , kt /)
       icount4 = (/ lx, ly, lz, 1  /)
       CALL disp_err(NF_GET_VARA_REAL(id_fil, id_var, istart4, icount4, X), &
            &        croutnm, cfil, cvar)
       !!
    END IF
    !!
    !! Closing file :
    IF (( kt == lt ).or.( kt == 0 ))  THEN
       CALL disp_err(NF_CLOSE(id_fil), croutnm, cfil, cvar)
       id_fil = 0 
       id_var = 0
    END IF
    !!
  END SUBROUTINE GETVAR_3D
  !!
  !!===============================================================================
  !!
  !!
  !!
  !!
  SUBROUTINE GETMASK(cfil, cvar, lx, ly, jlev, X)
    !!
    !!-----------------------------------------------------------------------
    !!  Get mask (variable 'cvar') from a netcdf file.
    !! - mask is stored in integer array X
    !!
    !! INPUT :
    !! -------
    !!          * cfil      : name of the input file          (character)
    !!          * cvar      : name of mask variable           (character)
    !!          * lx        : x dimension                       (integer)   
    !!          * ly        : y dimension                       (integer)
    !!          * jlev      : level to get (0 if no levels)     (integer)
    !!
    !! OUTPUT :
    !! --------
    !!          * X        :  2D array contening mask       (integer)
    !!
    !!
    !!  Author :            Laurent Brodeau, 2004
    !!  --------
    !!------------------------------------------------------------------------
    !!
    !!  
    CHARACTER(len=*),         INTENT(in)  :: cfil, cvar
    INTEGER,                    INTENT(in)  :: lx, ly, jlev
    INTEGER, DIMENSION(lx, ly), INTENT(out) :: X
    !!
    INTEGER :: ni, nj, nk, nt, icz
    !!
    croutnm = 'GETMASK'
    !!
    !!
    icz = 1
    IF ( jlev > 0 ) icz = jlev
    !!
    !!    -------------------------
    !!    Opening MASK netcdf file 
    !!    -------------------------  
    CALL disp_err(NF_OPEN(cfil, NF_NOWRITE, nfid), croutnm, cfil, cvar)
    !!
    ! Querrying dimensions on variable :
    CALL disp_err(NF_INQ_VARID(nfid, cvar, nvid), croutnm, cfil, cvar)
    !!
    CALL disp_err(NF_INQ_VARNDIMS(nfid,nvid,n_dims), croutnm, cfil, cvar)
    !!
    !!
    !!
    !!
    IF ( n_dims == 4 ) THEN
       !!
       CALL disp_err(NF_INQ_VARDIMID(nfid, nvid, id_dims4), croutnm, cfil, cvar)
       CALL disp_err(NF_INQ_DIMLEN(nfid, id_dims4(1), ni), croutnm, cfil, cvar)  
       CALL disp_err(NF_INQ_DIMLEN(nfid, id_dims4(2), nj), croutnm, cfil, cvar)
       CALL disp_err(NF_INQ_DIMLEN(nfid, id_dims4(3), nk), croutnm, cfil, cvar)
       CALL disp_err(NF_INQ_DIMLEN(nfid, id_dims4(4), nt), croutnm, cfil, cvar)
       !!
       IF ( (ni /= lx) .OR. (nj /= ly) ) THEN
          PRINT *, 'ERROR : data and mask file dont have same dimensions!'
          STOP
       ENDIF
       !!   
       !!  Filling mask(:,:) array with data values :   
       istart4 = (/ 1  , 1  , icz , 1 /)
       icount4 = (/ ni , nj , 1 , 1 /)
       !!
       CALL disp_err(NF_GET_VARA_INT(nfid,nvid,istart4,icount4,X), &
            &        croutnm, cfil, cvar)
       !!
       CALL disp_err(NF_CLOSE(nfid), croutnm, cfil, cvar)
       !!
       !!
       !!
    ELSE IF ( n_dims == 3 ) THEN
       !!***********************   
       !!
       CALL disp_err(NF_INQ_VARDIMID(nfid, nvid, id_dims3), croutnm, cfil, cvar)
       CALL disp_err(NF_INQ_DIMLEN  (nfid, id_dims3(1), ni),  croutnm, cfil, cvar)  
       CALL disp_err(NF_INQ_DIMLEN  (nfid, id_dims3(2), nj),  croutnm, cfil, cvar)
       CALL disp_err(NF_INQ_DIMLEN  (nfid, id_dims3(3), nt),  croutnm, cfil, cvar)
       !!
       IF ( (ni /= lx) .OR. (nj /= ly) ) THEN
          PRINT *, 'ERROR : data and mask file dont have same dimensions!'
          STOP
       ENDIF
       !!   
       !!  Filling mask(:,:) array with data values :   
       istart3 = (/ 1  , 1  , icz /)
       icount3 = (/ ni , nj , 1   /)
       !!
       CALL disp_err(NF_GET_VARA_INT(nfid,nvid,istart3,icount3,X), &
            &        croutnm, cfil, cvar)
       !!
       CALL disp_err(NF_CLOSE(nfid), croutnm, cfil, cvar)
       !!
       !!
       !!
    ELSE IF ( n_dims == 2 ) THEN
       !!
       IF ( jlev /= 0 ) THEN
          PRINT *, 'BAD! You want to get the land-sea mask at level', jlev
          PRINT *, 'BUT : you mask is 2D!'
          STOP
       END IF
       !!
       CALL disp_err(NF_INQ_VARDIMID(nfid,nvid,id_dims2), croutnm, cfil, cvar)
       CALL disp_err(NF_INQ_DIMLEN(nfid,id_dims2(1),ni), croutnm, cfil, cvar)  
       CALL disp_err(NF_INQ_DIMLEN(nfid,id_dims2(2),nj), croutnm, cfil, cvar)
       !!
       IF ( (ni /= lx) .OR. (nj /= ly) ) THEN
          PRINT *, 'ERROR : data and mask file dont have same spacial dimensions!'
          STOP
       ENDIF
       !!
       CALL disp_err(NF_GET_VAR_INT(nfid,nvid,X), croutnm, cfil, cvar)
       !!
       CALL disp_err(NF_CLOSE(nfid), croutnm, cfil, cvar)
       !!
       !!
       !!
    ELSE IF ( (n_dims /= 2).and.(n_dims /= 3).and.(n_dims /= 4) ) THEN
       PRINT *,'Error : mask should have 2, 3 or 4 dimension!'
       STOP
    ENDIF
    !!
  END SUBROUTINE GETMASK
  !!
  !!===============================================================================
  !!
  !!
  !!
  SUBROUTINE GETMASK_3D(cfil, cvar, lx, ly, lz, X)
    !!
    !!-----------------------------------------------------------------------
    !!  Get mask (variable 'cvar') from a netcdf file.
    !! - mask is stored in integer array X
    !!
    !! INPUT :
    !! -------
    !!          * cfil      : name of the input file          (character)
    !!          * cvar      : name of mask variable           (character)
    !!          * lx        : x dimension                       (integer)   
    !!          * ly        : y dimension                       (integer)
    !!          * lz        : z dimension                       (integer)
    !!
    !! OUTPUT :
    !! --------
    !!          * X        :  3D array contening mask       (integer)
    !!
    !!
    !!  Author :            Laurent Brodeau, 2004
    !!  --------
    !!------------------------------------------------------------------------
    !!
    !!  
    CHARACTER(len=*),         INTENT(in)  :: cfil, cvar
    INTEGER,                    INTENT(in)  :: lx, ly, lz
    INTEGER, DIMENSION(lx, ly, lz), INTENT(out) :: X
    !!
    INTEGER :: ni, nj, nk, nt
    !!
    croutnm = 'GETMASK_3D'
    !!
    !!    -------------------------
    !!    Opening MASK netcdf file 
    !!    -------------------------  
    CALL disp_err(NF_OPEN(cfil, NF_NOWRITE, nfid), croutnm, cfil, cvar)
    !!
    ! Querrying dimensions on variable :
    CALL disp_err(NF_INQ_VARID(nfid, cvar, nvid), croutnm, cfil, cvar)
    !!
    CALL disp_err(NF_INQ_VARNDIMS(nfid,nvid,n_dims), croutnm, cfil, cvar)
    !!
    !!
    !!
    !!
    IF ( n_dims == 4 ) THEN
       !!
       CALL disp_err(NF_INQ_VARDIMID(nfid, nvid, id_dims4), croutnm, cfil, cvar)
       CALL disp_err(NF_INQ_DIMLEN(nfid, id_dims4(1), ni), croutnm, cfil, cvar)  
       CALL disp_err(NF_INQ_DIMLEN(nfid, id_dims4(2), nj), croutnm, cfil, cvar)
       CALL disp_err(NF_INQ_DIMLEN(nfid, id_dims4(3), nk), croutnm, cfil, cvar)
       CALL disp_err(NF_INQ_DIMLEN(nfid, id_dims4(4), nt), croutnm, cfil, cvar)
       !!
       IF ( (ni /= lx) .OR. (nj /= ly) .OR. (nk /= lz) ) THEN
          PRINT *, 'ERROR : data and mask file dont have same dimensions!'
          STOP
       ENDIF
       !!   
       !!  Filling mask(:,:) array with data values :   
       istart4 = (/ 1  , 1  , 1 , 1 /)
       icount4 = (/ ni , nj , nk , 1 /)
       !!
       CALL disp_err(NF_GET_VARA_INT(nfid, nvid, istart4, icount4, X), croutnm, cfil, cvar)
       !!
       CALL disp_err(NF_CLOSE(nfid), croutnm, cfil, cvar)
       !!
       !!
       !!
       !!
    ELSE IF ( n_dims == 3 ) THEN
       !!***********************   
       !!
       CALL disp_err(NF_INQ_VARDIMID(nfid, nvid, id_dims3), croutnm, cfil, cvar)
       CALL disp_err(NF_INQ_DIMLEN  (nfid, id_dims3(1), ni),  croutnm, cfil, cvar)  
       CALL disp_err(NF_INQ_DIMLEN  (nfid, id_dims3(2), nj),  croutnm, cfil, cvar)
       CALL disp_err(NF_INQ_DIMLEN  (nfid, id_dims3(3), nt),  croutnm, cfil, cvar)
       !!
       IF ( (ni /= lx) .OR. (nj /= ly) .OR. (nk /= lz) ) THEN
          PRINT *, 'ERROR : data and mask file dont have same dimensions!'
          STOP
       ENDIF
       !!   
       !!  Filling mask(:,:) array with data values :   
       istart3 = (/ 1  , 1  , 1   /)
       icount3 = (/ ni , nj , nk  /)
       !!
       CALL disp_err(NF_GET_VARA_INT(nfid, nvid, istart3, icount3, X), croutnm, cfil, cvar)
       !!
       CALL disp_err(NF_CLOSE(nfid), croutnm, cfil, cvar)
       !!
       !!
    ELSE IF ( (n_dims /= 3).and.(n_dims /= 4) ) THEN
       PRINT *,'Error : mask should have 2, 3 or 4 dimension!'
       STOP
    ENDIF
    !!
  END SUBROUTINE GETMASK_3D
  !!
  !!===============================================================================
  !!
  !!
  !!
  !!
  SUBROUTINE P2D_T(id_fil, id_var, lx, ly, lt, lct, vlon, vlat, vtime, &
       &         x2d, cfil, cvarlon, cvarlat, cvartime, cvar, cunit, &
       &         cln, vflag, cun_t, cvlonout, cvlatout)
    !!
    !!-----------------------------------------------------------------------------
    !! This routine prints a 3D array to a netcdf file with regular
    !! coordinates as 2 new variables 1D : 'lon(lon)' and 'lat(lat)'
    !!
    !! INPUT :
    !! -------
    !!        id_fil = ID of the file (takes its value on the first call)
    !!        id_var = ID of the variable //
    !!        lx    = x dimension of array to plot              [integer]
    !!        ly    = y dimension of array to plot              [integer]
    !!        lt    = t dimension of array to plot              [integer]
    !!        lct   = current time step                         [integer]
    !!        vlon  = 1D array (lx) of longitudes               [real]
    !!        vlat  = 1D array (ly) of latitudes                [real]
    !!        vtime = 1D array (lt) of times                    [real]
    !!        x2d   = 2D snap of 3D array at time jt to write   [real]
    !!        cfil  = name of the output file                   [character]
    !!        cvarlon = name of the longitude                   [character]
    !!        cvarlat = name of the latitude                    [character]
    !!        cvartime = name of time                           [character]
    !!        cvar  = name of the variable                      [character]
    !!        cunit = unit of the output variable               [character]
    !!        cln   = long name of variable to be treated       [character]
    !!        vflag = flag value                                [real]
    !!
    !!        cun_t = unit for time                         |OPTIONAL|  [character]
    !!        cvlonout = name for lon array in output file  |OPTIONAL|  [character]
    !!        cvlatout = name for lat array in output file  |OPTIONAL|  [character]
    !!------------------------------------------------------------------------------
    !!
    INTEGER, INTENT(inout) :: id_fil, id_var
    INTEGER, INTENT(in)    :: lx, ly, lt, lct
    !!
    REAL(wpez), DIMENSION(lx),       INTENT(in) :: vlon
    REAL(wpez), DIMENSION(ly),       INTENT(in) :: vlat
    REAL(wpez), DIMENSION(lt),       INTENT(in) :: vtime
    REAL(wpez), DIMENSION(lx,ly),    INTENT(in) :: x2d
    !!
    CHARACTER(len=*), INTENT(in) :: cfil, cvarlon, cvarlat, cvartime, cvar, cunit, cln
    CHARACTER(len=*),  OPTIONAL  :: cun_t, cvlonout, cvlatout
    !!
    REAL(wpez),             INTENT(in) :: vflag
    !!
    !!
    croutnm = 'P2D_T'
    !!
    cvxout = trim(cvarlon) ; cvyout = trim(cvarlat)
    IF ( present(cvlonout) ) cvxout = trim(cvlonout)
    IF ( present(cvlatout) ) cvyout = trim(cvlatout)
    !!
    !!
    IF ( lct == 1 ) THEN
       !!
       !!      CREATE NETCDF OUTPUT FILE :
       CALL disp_err(NF_CREATE(cfil,NF_CLOBBER,id_fil), croutnm, cfil, cvar)
       !!
       !!      Create longitude dimension :
       CALL disp_err(NF_DEF_DIM(id_fil, trim(cvarlon), lx, id_x), croutnm, cfil, cvar)
       !!
       !       Create latitude dimension :
       CALL disp_err(NF_DEF_DIM(id_fil, trim(cvarlat), ly, id_y), croutnm, cfil, cvar)
       !!
       !!      Create record dimension :
       CALL disp_err(NF_DEF_DIM(id_fil, trim(cvartime), NF_UNLIMITED, id_t), &
            &        croutnm, cfil, cvar)
       !!
       id_dims3 = (/ id_x, id_y, id_t /)
       !!
       !!
       !!
       !!           VARIABLE TO PLOT
       !!           ----------------
       !! Longitude
       CALL disp_err(NF_DEF_VAR(id_fil, trim(cvxout),NF_REAL, 1, id_x, id_lon), &
            &        croutnm, cfil, cvar)
       !!
       !! Latitude
       CALL disp_err(NF_DEF_VAR(id_fil, trim(cvyout),NF_REAL, 1, id_y, id_lat), &
            &        croutnm, cfil, cvar)
       !!
       !! Time
       CALL disp_err(NF_DEF_VAR(id_fil,trim(cvartime),NF_REAL, 1, id_t, id_time), &
            &        croutnm, cfil, cvar)
       !!
       !! Variable
       CALL disp_err(NF_DEF_VAR(id_fil,trim(cvar),NF_REAL,3,id_dims3,id_var),     &
            &      croutnm, cfil, cvar)
       !!
       !!
       !!
       !!          ATTRIBUTES
       !!          ----------
       !!
       !!      For longitude :
       CALL disp_err(NF_PUT_ATT_TEXT(id_fil, id_lon, 'units', 12,'degrees_east'), &
            &      croutnm, cfil, cvar)
       !!
       rrange = (/ minval(vlon) , maxval(vlon) /)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_lon,'valid_min',NF_REAL,1,rrange(1)), &
            &        croutnm, cfil, cvar)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_lon,'valid_max',NF_REAL,1,rrange(2)), &
            &        croutnm, cfil, cvar)
       !!
       !!      For latitude :
       CALL disp_err(NF_PUT_ATT_TEXT(id_fil, id_lat, 'units', 13,'degrees_north'),  &
            &        croutnm, cfil, cvar) 
       !!
       rrange = (/ minval(vlat) , maxval(vlat) /)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_lat,'valid_min',NF_REAL,1,rrange(1)), &
            &        croutnm, cfil, cvar)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_lat,'valid_max',NF_REAL,1,rrange(2)), &
            &        croutnm, cfil, cvar)
       !!
       !!      For time
       cu = 'unknown'
       IF ( present(cun_t) ) cu = cun_t
       ils = len_trim(trim(cu))
       CALL disp_err(NF_PUT_ATT_TEXT(id_fil, id_time, 'units', ils, trim(cu)), &
            &        croutnm, cfil, cvar)
       rrange = (/ minval(vtime), maxval(vtime) /)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_time,'valid_min',NF_REAL,1,rrange(1)), &
            &        croutnm, cfil, cvar)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_time,'valid_max',NF_REAL,1,rrange(2)), &
            &        croutnm, cfil, cvar)
       !!
       !!
       !! FOR VARIABLE :
       !!
       !! Long name
       ils=LEN_TRIM(trim(cln))
       CALL disp_err(NF_PUT_ATT_TEXT(id_fil, id_var, 'long_name', ils, trim(cln)), &
            &       croutnm, cfil, cvar)
       !!
       !! Units
       ils=LEN_TRIM(cunit) 
       CALL disp_err(NF_PUT_ATT_TEXT(id_fil, id_var, 'units', ils, trim(cunit) ),  &
            &       croutnm, cfil, cvar)     
       !!
       !! Missing value
       IF ( vflag /= 0. ) &
            & CALL disp_err(NF_PUT_ATT_REAL(id_fil, id_var, 'missing_value', &
            &              NF_REAL, 1, vflag), croutnm, cfil, cvar)
       !!
       !!
       !! Calculating min without the mask vflag :
       IF ( vflag /= 0.) THEN
          rmin = 1.E4
          DO ji=1, lx
             DO jj=1, ly
                IF ((x2d(ji,jj) <= rmin).and.(x2d(ji,jj) > vflag)) rmin = x2d(ji,jj)
             END DO
          END DO
       ELSE
          rmin = minval(x2d)
       END IF
       !!
       rrange = (/ rmin , maxval(x2d) /)
       !!
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_var,'valid_range',NF_REAL,2,rrange), &
            &        croutnm, cfil, cvar)
       !!
       !! Global attributes
       ils = LEN_TRIM(cabout)
       CALL disp_err(NF_PUT_ATT_TEXT(id_fil, NF_GLOBAL, 'About', ils, trim(cabout)), &
            &        croutnm, cfil, cvar)
       !!
       !!
       !!           END OF DEFINITION
       !!           -----------------
       CALL disp_err(NF_ENDDEF(id_fil), croutnm, cfil, cvar)
       !!
       !!       Write longitude variable :
       CALL disp_err(NF_PUT_VAR_REAL(id_fil, id_lon, vlon), croutnm, cfil, cvar)
       !!
       !!       Write latitude variable :
       CALL disp_err(NF_PUT_VAR_REAL(id_fil, id_lat, vlat), croutnm, cfil, cvar)
       !!
    END IF
    !!
    !!          WRITE TIME
    CALL disp_err(NF_PUT_VARA_REAL(id_fil, id_t, lct, 1, vtime(lct)), &
         &        croutnm, cfil, cvar)
    !!
    !!          WRITE VARIABLE
    istart3 = (/ 1, 1, lct /)
    icount3 = (/ lx, ly, 1 /)
    CALL disp_err(NF_PUT_VARA_REAL(id_fil, id_var, istart3, icount3, x2d ), &
         &        croutnm, cfil, cvar)
    !!
    IF ( lct == lt ) CALL disp_err(NF_CLOSE(id_fil), croutnm, cfil, cvar)
    !!
    !!
  END SUBROUTINE P2D_T
  !!
  !!
  !!
  !!
  !!
  SUBROUTINE P2D_T_IRR(id_fil, id_var, lx, ly, lt, lct, rlon, rlat, vtime, x2d,     &
       &             cfil, cvarlon, cvarlat, cvartime, cvar, cunit, cln, vflag, cun_t)
    !!
    !!
    !! INPUT :
    !! -------
    !!        id_fil = ID of the file (takes its value on the first call)
    !!        id_var = ID of the variable //
    !!        lx    = x dimension of array to plot             [integer]
    !!        ly    = y dimension of array to plot             [integer]
    !!        lt    = t dimension of array to plot             [integer]
    !!        lct   = current time step                         [integer]
    !!        rlon  = 2D array of longitude                     [real]
    !!        rlat  = 2D array of latitude                      [real]
    !!        vtime  = time array                               [array 1D]
    !!        x2d   = 2D snap of 3D array at time jt to write   [real]
    !!        cfil  = name of the output file                   [character]
    !!        cvarlon = name of longitude                       [character]
    !!        cvarlat = name of latitude                        [character]
    !!        cvartime = name of time                           [character]
    !!        cvar  = name of the variable                      [character]
    !!        vflag = flag value or "0."                        [real]
    !!        cunit  = unit for treated variable                [character]
    !!        cln = long-name for treated variable        [character]
    !!        cfilename = actual file name (without path!)      [character]
    !!
    !!        cun_t = unit for time                 |OPTIONAL|  [character]
    !!
    !!--------------------------------------------------------------------------
    !!
    !!
    INTEGER, INTENT(inout) :: id_fil, id_var
    INTEGER, INTENT(in)    :: lx, ly, lt, lct
    !!
    REAL(wpez), DIMENSION(lx,ly), INTENT(in) :: rlat, rlon, x2d
    REAL(wpez), DIMENSION(lt),    INTENT(in) :: vtime
    !!
    CHARACTER(len=*), INTENT(in) :: cfil, cvarlon, cvarlat, cvartime, cvar, cunit, cln
    CHARACTER(len=*),  OPTIONAL, INTENT(in) :: cun_t
    !!
    REAL(wpez),         INTENT(in) :: vflag
    !!
    !!
    croutnm = 'P2D_T_IRR'
    !!
    !!
    IF ( lct == 1 ) THEN
       !!
       !!
       !! Opening mesh file for grid quest :
       !! ----------------------------------
       !!
       !!           CREATE NETCDF OUTPUT FILE :
       CALL disp_err(NF_CREATE(cfil, NF_CLOBBER, id_fil), croutnm, cfil, cvar)
       !!
       !!             DIMMENSIONS 
       !!       Create longitude dimension :
       CALL disp_err(NF_DEF_DIM(id_fil, 'x', lx, id_x), croutnm, cfil, cvar)
       !!
       !       Create latitude dimension :
       CALL disp_err(NF_DEF_DIM(id_fil, 'y', ly, id_y), croutnm, cfil, cvar)
       !!
       !!      Create time record dimension :
       CALL disp_err(NF_DEF_DIM(id_fil, trim(cvartime), NF_UNLIMITED, id_t), &
            &        croutnm, cfil, cvar)
       !CALL disp_err(NF_DEF_DIM(id_fil, trim(cvartime), lt, id_t), &
       !    &        croutnm, cfil, cvar)
       !!
       !!
       id_dims2 = (/ id_x , id_y /)
       id_dims3 = (/ id_x , id_y , id_t /)
       !!
       !!
       !!           VARIABLES TO PLOT
       !!       Create longitude variable :
       CALL disp_err(NF_DEF_VAR(id_fil,trim(cvarlon),NF_REAL,2,id_dims2,id_lon), &
            &      croutnm, cfil, cvar)
       !!
       !!       Create latitude variable :
       CALL disp_err(NF_DEF_VAR(id_fil,trim(cvarlat),NF_REAL,2,id_dims2,id_lat), &
            &       croutnm, cfil, cvar)  
       !!
       !! Time
       CALL disp_err(NF_DEF_VAR(id_fil,trim(cvartime),NF_REAL,1,id_t,id_time), &
            &        croutnm, cfil, cvar)
       !!
       !! Variable
       CALL disp_err(NF_DEF_VAR(id_fil,trim(cvar),NF_REAL,3,id_dims3,id_var), &
            &       croutnm, cfil, cvar)
       !!
       !!
       !!          ATTRIBUTES
       !!
       !!      For longitude :
       CALL disp_err(NF_PUT_ATT_TEXT(id_fil, id_lon, 'units', 12,'degrees_east'), &
            &      croutnm, cfil, cvar)
       !!
       rrange = (/ minval(rlon) , maxval(rlon) /)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_lon,'valid_min',NF_REAL,1,rrange(1)) &
            &       , croutnm, cfil, cvar)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_lon,'valid_max',NF_REAL,1,rrange(2)) &
            &       , croutnm, cfil, cvar)
       !!
       !!      For latitude :
       CALL disp_err(NF_PUT_ATT_TEXT(id_fil, id_lat, 'units', 13,'degrees_north')  &
            &       , croutnm, cfil, cvar) 
       !!
       rrange = (/ minval(rlat) , maxval(rlat) /)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_lat,'valid_min',NF_REAL,1,rrange(1)) &
            &       , croutnm, cfil, cvar)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_lat,'valid_max',NF_REAL,1,rrange(2)) &
            &       , croutnm, cfil, cvar)
       !!
       !!      For time
       cu = 'unknown'
       IF ( present(cun_t) ) cu = cun_t
       ils = len_trim(trim(cu))
       CALL disp_err(NF_PUT_ATT_TEXT(id_fil, id_time, 'units', ils, trim(cu)), &
            &        croutnm, cfil, cvar)
       rrange = (/ minval(vtime), maxval(vtime) /)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_time,'valid_min',NF_REAL,1,rrange(1)),&
            &        croutnm, cfil, cvar)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_time,'valid_max',NF_REAL,1,rrange(2)),&
            &          croutnm, cfil, cvar)
       !!
       !!
       !!  VARIABLE ATTRIBUTES
       !!
       !! Long name
       ils=LEN_TRIM(trim(cln))
       CALL disp_err(NF_PUT_ATT_TEXT(id_fil, id_var, 'long_name', ils, trim(cln)), &
            &       croutnm, cfil, cvar)
       !!
       !! Units
       ils = LEN_TRIM(cunit) 
       CALL disp_err(NF_PUT_ATT_TEXT(id_fil, id_var, 'units', ils, trim(cunit) ),  &
            &       croutnm, cfil, cvar)     
       !!
       !! Missing value
       IF ( vflag /= 0.) &
            & CALL disp_err(NF_PUT_ATT_REAL(id_fil, id_var, 'missing_value', &
            &              NF_REAL, 1, vflag), croutnm, cfil, cvar)
       !!
       !! Calculating min without the mask vflag :
       IF ( vflag /= 0.) THEN
          rmin = 1.E4
          DO ji=1, lx
             DO jj=1, ly
                IF ((x2d(ji,jj) <= rmin).and.(x2d(ji,jj) > vflag)) rmin = x2d(ji,jj)
             END DO
          END DO
       ELSE
          rmin = minval(x2d)
       END IF
       !!
       rrange = (/ rmin , maxval(x2d) /)
       !!
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_var,'valid_range',NF_REAL,2,rrange), &
            &        croutnm, cfil, cvar)
       !!
       !!
       !! Global attributes
       ils = LEN_TRIM(cabout)
       CALL disp_err(NF_PUT_ATT_TEXT(id_fil, NF_GLOBAL, 'About', ils, trim(cabout))  &
            &       , croutnm, cfil, cvar)
       !!
       !!
       !!           END OF DEFINITION
       !!           -----------------
       CALL disp_err(NF_ENDDEF(id_fil), croutnm, cfil, cvar)
       !!
       !!
       !!          WRITE COORDINATES
       !!          -----------------
       !!       Write longitude variable :
       CALL disp_err(NF_PUT_VAR_REAL(id_fil, id_lon, rlon), croutnm, cfil, cvar)
       !!
       !!       Write latitude variable :
       CALL disp_err(NF_PUT_VAR_REAL(id_fil, id_lat, rlat), croutnm, cfil, cvar)
       !!
       !!       Write time variable :
       CALL disp_err(NF_PUT_VARA_REAL(id_fil, id_time, 1, lt, vtime), &
            &        croutnm, cfil, cvar)
       !!
    END IF
    !!
    !!
    !!          WRITE VARIABLE
    istart3 = (/ 1,  1, lct /)
    icount3 = (/ lx, ly, 1 /)
    CALL disp_err(NF_PUT_VARA_REAL(id_fil, id_var, istart3, icount3, x2d), &
         &        croutnm, cfil, cvar)
    !!
    IF ( lct == lt ) CALL disp_err(NF_CLOSE(id_fil), croutnm, cfil, cvar)
    !!
    !!
  END SUBROUTINE P2D_T_IRR
  !!
  !!
  !!
  !!
  SUBROUTINE P3D_T(id_fil, id_var, lx, ly, lz, lt, lct, vlon, vlat, vdpth, vtime, &
       &           x3d, cfil, cvarlon, cvarlat, cvardpth, cvartime, cvar, cunit, &
       &           cln, vflag, cun_z, cun_t)
    !!
    !!
    !! INPUT :
    !! -------
    !!        id_fil = ID of the file (takes its value on the first call)
    !!        id_var = ID of the variable //
    !!        lx    = x dimension of array to plot             [integer]
    !!        ly    = y dimension of array to plot             [integer]
    !!        lz    = z dimension of array to plot             [integer]
    !!        lt    = t dimension of array to plot             [integer]
    !!        lct   = current time step                         [integer]
    !!        vlon  = 1D array of longitude                     [real]
    !!        vlat  = 1D array of latitude                      [real]
    !!        vdpth = depth array                               [array 1D]
    !!        vtime = time array                                [array 1D]
    !!        x3d   = 3D snapshot at time jt to write           [real]
    !!        cfil  = name of the output file                   [character]
    !!        cvarlon = name of longitude                       [character]
    !!        cvarlat = name of latitude                        [character]
    !!        cvardpth = name of depth                          [character]
    !!        cvartime = name of time                           [character]
    !!        cvar  = name of the variable                      [character]
    !!        vflag = flag value or "0."                        [real]
    !!        cunit  = unit for treated variable                [character]
    !!        cln = long-name for treated variable              [character]
    !!        cfilename = actual file name (without path!)      [character]
    !!
    !!        cun_z = unit for depth                |OPTIONAL|  [character]
    !!        cun_t = unit for time                 |OPTIONAL|  [character]
    !!--------------------------------------------------------------------------
    !!
    !!
    INTEGER, INTENT(inout) :: id_fil, id_var
    INTEGER, INTENT(in)    :: lx, ly, lz, lt, lct
    !!
    REAL(wpez), DIMENSION(lx), INTENT(in) :: vlon
    REAL(wpez), DIMENSION(ly), INTENT(in) :: vlat
    REAL(wpez), DIMENSION(lx,ly,lz), INTENT(in) :: x3d
    REAL(wpez), DIMENSION(lz),    INTENT(in) :: vdpth
    REAL(wpez), DIMENSION(lt),    INTENT(in) :: vtime
    !!
    CHARACTER(len=*), INTENT(in) :: cfil, cvarlon, cvarlat, cvardpth, cvartime, cvar, cunit, cln
    CHARACTER(len=*),  OPTIONAL, INTENT(in) :: cun_z, cun_t
    !!
    REAL(wpez),         INTENT(in) :: vflag
    !!
    !!
    croutnm = 'P3D_T'
    !!
    !!
    IF ( lct == 1 ) THEN
       !!
       !!
       !! Opening mesh file for grid quest :
       !! ----------------------------------
       !!
       !!           CREATE NETCDF OUTPUT FILE :
       CALL disp_err(NF_CREATE(cfil, NF_CLOBBER, id_fil), croutnm, cfil, cvar)
       !!
       !!             DIMMENSIONS 
       !!       Create longitude dimension :
       CALL disp_err(NF_DEF_DIM(id_fil, cvarlon, lx, id_x), croutnm, cfil, cvar)
       !!
       !!       Create latitude dimension :
       CALL disp_err(NF_DEF_DIM(id_fil, cvarlat, ly, id_y), croutnm, cfil, cvar)
       !!
       !!       Create depth dimension :
       CALL disp_err(NF_DEF_DIM(id_fil, cvardpth, lz, id_z), croutnm, cfil, cvar)
       !!
       !!      Create record dimension :
       CALL disp_err(NF_DEF_DIM(id_fil, trim(cvartime), NF_UNLIMITED, id_t), &
            &        croutnm, cfil, cvar)
       !!
       id_dims4 = (/ id_x , id_y , id_z, id_t /)
       !!
       !!           VARIABLES TO PLOT
       !!       Create longitude variable :
       CALL disp_err(NF_DEF_VAR(id_fil, trim(cvarlon), NF_REAL, 1, id_x, id_lon),  &
          &      croutnm, cfil, cvar)
       !!
       !!       Create latitude variable :
       CALL disp_err(NF_DEF_VAR(id_fil, trim(cvarlat), NF_REAL, 1, id_y, id_lat), &
            &      croutnm, cfil, cvar)  
       !!
       !!       Create depth variable :
       CALL disp_err(NF_DEF_VAR(id_fil, trim(cvardpth), NF_REAL, 1, id_z, id_dpth), &
            &      croutnm, cfil, cvar)  
       !!
       !! Time
       CALL disp_err(NF_DEF_VAR(id_fil, trim(cvartime),  NF_REAL, 1, id_t, id_time), &
            &      croutnm, cfil, cvar)
       !!
       !! Variable
       CALL disp_err(NF_DEF_VAR(id_fil, trim(cvar), NF_REAL, 4, id_dims4, id_var), &
            &      croutnm, cfil, cvar)
       !!
       !!          ATTRIBUTES
       !!
       !!      For longitude :
       CALL disp_err(NF_PUT_ATT_TEXT(id_fil, id_lon, 'units', 12,'degrees_east'), &
            &      croutnm, cfil, cvar)
       !!
       rrange = (/ minval(vlon) , maxval(vlon) /)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_lon,'valid_min',NF_REAL,1,rrange(1)), &
            &        croutnm, cfil, cvar)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_lon,'valid_max',NF_REAL,1,rrange(2)), &
            &        croutnm, cfil, cvar)
       !!
       !!      For latitude :
       CALL disp_err(NF_PUT_ATT_TEXT(id_fil, id_lat, 'units', 13,'degrees_north'),  &
            &        croutnm, cfil, cvar) 
       !!
       rrange = (/ minval(vlat) , maxval(vlat) /)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_lat,'valid_min',NF_REAL,1,rrange(1)), &
            &        croutnm, cfil, cvar)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_lat,'valid_max',NF_REAL,1,rrange(2)), &
            &        croutnm, cfil, cvar)
       !!
       !!      For depth :
       cu = 'unknown'
       IF ( present(cun_z) ) cu = cun_z
       ils = len_trim(trim(cu))
       CALL disp_err(NF_PUT_ATT_TEXT(id_fil, id_dpth, 'units', ils, trim(cu)), &
            &        croutnm, cfil, cvar) 
       !!
       rrange = (/ minval(vdpth) , maxval(vdpth) /)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_dpth,'valid_min',NF_REAL,1,rrange(1)), &
            &       croutnm, cfil, cvar)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_dpth,'valid_max',NF_REAL,1,rrange(2)), &
            &       croutnm, cfil, cvar)
       !!
       !!      For time
       cu = 'unknown'
       IF ( present(cun_t) ) cu = cun_t
       ils = len_trim(trim(cu))
       CALL disp_err(NF_PUT_ATT_TEXT(id_fil, id_time, 'units', ils, trim(cu)), &
            &        croutnm, cfil, cvar)
       rrange = (/ minval(vtime), maxval(vtime) /)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_time,'valid_min',NF_REAL,1,rrange(1)), &
            &        croutnm, cfil, cvar)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_time,'valid_max',NF_REAL,1,rrange(2)), &
            &        croutnm, cfil, cvar)
       !!
       !!
       !!  VARIABLE ATTRIBUTES
       !!
       !! Long name
       ils=LEN_TRIM(trim(cln))
       CALL disp_err(NF_PUT_ATT_TEXT(id_fil, id_var, 'long_name', ils, trim(cln)), &
            &       croutnm, cfil, cvar)
       !!
       !! Units
       ils = LEN_TRIM(cunit) 
       CALL disp_err(NF_PUT_ATT_TEXT(id_fil, id_var, 'units', ils, trim(cunit) ),  &
            &       croutnm, cfil, cvar)     
       !!
       !! Missing value
       IF ( vflag /= 0.) &
            & CALL disp_err(NF_PUT_ATT_REAL(id_fil, id_var, 'missing_value', &
            &              NF_REAL, 1, vflag), croutnm, cfil, cvar)
       !!
       !! Calculating min without the mask vflag :
       IF ( vflag /= 0.) THEN
          rmin = 1.E4
          DO ji=1, lx
             DO jj=1, ly
                IF ((x3d(ji,jj,1) <= rmin).and.(x3d(ji,jj,1) > vflag)) &
                     rmin = x3d(ji,jj,1)
             END DO
          END DO
       ELSE
          rmin = minval(x3d)
       END IF
       !!
       rrange = (/ rmin , maxval(x3d) /)
       !!
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_var,'valid_range',NF_REAL,2,rrange), &
            &        croutnm, cfil, cvar)
       !!
       !!
       !!
       !! Global attributes
       ils = LEN_TRIM(cabout)
       CALL disp_err(NF_PUT_ATT_TEXT(id_fil, NF_GLOBAL, 'About', ils, trim(cabout)), &
            &        croutnm, cfil, cvar)
       !!
       !!
       !!           END OF DEFINITION
       !!           -----------------
       CALL disp_err(NF_ENDDEF(id_fil), croutnm, cfil, cvar)
       !!
       !!
       !!          WRITE COORDINATES
       !!          -----------------
       !!       Write longitude variable :
       CALL disp_err(NF_PUT_VAR_REAL(id_fil, id_lon, vlon), croutnm, cfil, cvar)
       !!
       !!       Write latitude variable :
       CALL disp_err(NF_PUT_VAR_REAL(id_fil, id_lat, vlat), croutnm, cfil, cvar)
       !!
       !!       Write depth variable :
       CALL disp_err(NF_PUT_VAR_REAL(id_fil, id_dpth, vdpth), croutnm, cfil, cvar)
       !!
       !!       Write time variable :
       CALL disp_err(NF_PUT_VARA_REAL(id_fil, id_time, 1, lt, vtime), &
            &        croutnm, cfil, cvar)
       !!
    END IF
    !!
    !!
    !!          WRITE VARIABLE
    istart4 = (/ 1,  1, 1, lct /)
    icount4 = (/ lx, ly, lz, 1 /)
    CALL disp_err(NF_PUT_VARA_REAL(id_fil, id_var, istart4, icount4, x3d), &
         &        croutnm, cfil, cvar)
    !!
    !!
    IF ( lct == lt ) CALL disp_err(NF_CLOSE(id_fil), croutnm, cfil, cvar)
    !!
    !!
  END SUBROUTINE P3D_T
  !!
  !!
  !!
  !!
  SUBROUTINE P3D_T_IRR(id_fil, id_var, linit, lx, ly, lz, lt, lct, rlon, rlat, &
       &               vdpth, &
       &               vtime, x3d, cfil, cvarlon, cvarlat, cvardpth, cvartime, &
       &               cvar, cunit, cln, vflag, cun_z, cun_t)
    !!
    !! INPUT :
    !! -------
    !!        id_fil = ID of the file (takes its value on the first call)
    !!        id_var = ID of the variable //
    !!        linit = flag for initializing file (0/1; no/yes)
    !!        lx    = x dimension of array to plot             [integer]
    !!        ly    = y dimension of array to plot             [integer]
    !!        lz    = z dimension of array to plot             [integer]
    !!        lt    = t dimension of array to plot             [integer]
    !!        lct   = current time step                         [integer]
    !!        rlon  = 2D array of longitude                     [real]
    !!        rlat  = 2D array of latitude                      [real]
    !!        vdpth = depth array                               [array 1D]
    !!        vtime = time array                                [array 1D]
    !!        x3d   = 3D snapshot at time jt to write           [real]
    !!        cfil  = name of the output file                   [character]
    !!        cvarlon = name of longitude                       [character]
    !!        cvarlat = name of latitude                        [character]
    !!        cvardpth = name of depth                          [character]
    !!        cvartime = name of time                           [character]
    !!        cvar  = name of the variable                      [character]
    !!        vflag = flag value or "0."                        [real]
    !!        cunit  = unit for treated variable                [character]
    !!        cln = long-name for treated variable              [character]
    !!        cfilename = actual file name (without path!)      [character]
    !!
    !!        cun_z = unit for depth                |OPTIONAL|  [character]
    !!        cun_t = unit for time                 |OPTIONAL|  [character]
    !!
    !!--------------------------------------------------------------------------
    !!
    !!
    INTEGER, INTENT(inout) :: id_fil, id_var
    INTEGER, INTENT(in)    :: lx, ly, lz, lt, lct, linit
    !!
    REAL(wpez), DIMENSION(lx,ly), INTENT(in) :: rlat, rlon
    REAL(wpez), DIMENSION(lx,ly,lz), INTENT(in) :: x3d
    REAL(wpez), DIMENSION(lz),    INTENT(in) :: vdpth
    REAL(wpez), DIMENSION(lt),    INTENT(in) :: vtime
    !!
    CHARACTER(len=*), INTENT(in) :: cfil, cvarlon, cvarlat, cvardpth, cvartime, cvar, cunit, cln
    CHARACTER(len=*),  OPTIONAL, INTENT(in) :: cun_z, cun_t
    !!
    REAL(wpez),         INTENT(in) :: vflag
    !!
    !!
    croutnm = 'P3D_T_IRR'
    !!
    !!
    IF ( lct == 1 ) THEN
      IF (linit == 1) THEN
       !!
       !!
       !! Opening mesh file for grid quest :
       !! ----------------------------------
       !!
       !!           CREATE NETCDF OUTPUT FILE :
       CALL disp_err(NF_CREATE(cfil, NF_CLOBBER, id_fil), croutnm, cfil, cvar)
       !!
       !!             DIMMENSIONS 
       !!       Create longitude dimension :
       CALL disp_err(NF_DEF_DIM(id_fil, 'x', lx, id_x), croutnm, cfil, cvar)
       !!
       !       Create latitude dimension :
       CALL disp_err(NF_DEF_DIM(id_fil, 'y', ly, id_y), croutnm, cfil, cvar)
       !!
       !       Create latitude dimension :
       CALL disp_err(NF_DEF_DIM(id_fil, trim(cvardpth), lz, id_z), croutnm, cfil, cvar)
       !!
       !!      Create record dimension :
       CALL disp_err(NF_DEF_DIM(id_fil, trim(cvartime), NF_UNLIMITED, id_t), &
            &        croutnm, cfil, cvar)
       !!
       !!
       id_dims2 = (/ id_x , id_y /)
       id_dims4 = (/ id_x , id_y , id_z, id_t /)
       !!
       !!
       !!           VARIABLES TO PLOT
       CALL disp_err(NF_DEF_VAR(id_fil,trim(cvarlon),NF_REAL,2,id_dims2,id_lon), &
            &        croutnm, cfil, cvar)
       !!
       CALL disp_err(NF_DEF_VAR(id_fil,trim(cvarlat),NF_REAL,2,id_dims2,id_lat),  &
            &        croutnm, cfil, cvar)  
       !!
       CALL disp_err(NF_DEF_VAR(id_fil,trim(cvardpth),NF_REAL,1,id_z,id_dpth),  &
            &        croutnm, cfil, cvar)  
       !!
       CALL disp_err(NF_DEF_VAR(id_fil,trim(cvartime), NF_REAL,1,id_t,id_time),  &
            &        croutnm, cfil, cvar)
       !!
       !!
       !!          ATTRIBUTES
       !!
       !!      For longitude :
       CALL disp_err(NF_PUT_ATT_TEXT(id_fil, id_lon, 'units', 12,'degrees_east'), &
            &      croutnm, cfil, cvar)
       !!
       rrange = (/ minval(rlon) , maxval(rlon) /)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_lon,'valid_min',NF_REAL,1,rrange(1)), &
            &        croutnm, cfil, cvar)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_lon,'valid_max',NF_REAL,1,rrange(2)), &
            &        croutnm, cfil, cvar)
       !!
       !!      For latitude :
       CALL disp_err(NF_PUT_ATT_TEXT(id_fil, id_lat, 'units', 13,'degrees_north'),  &
            &        croutnm, cfil, cvar) 
       !!
       rrange = (/ minval(rlat) , maxval(rlat) /)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_lat,'valid_min',NF_REAL,1,rrange(1)), &
            &        croutnm, cfil, cvar)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_lat,'valid_max',NF_REAL,1,rrange(2)), &
            &        croutnm, cfil, cvar)
       !!
       !!      For depth :
       cu = 'unknown'
       IF ( present(cun_z) ) cu = cun_z
       ils = len_trim(trim(cu))
       CALL disp_err(NF_PUT_ATT_TEXT(id_fil, id_dpth, 'units', ils, trim(cu)), &
            &        croutnm, cfil, cvar) 
       !!
       rrange = (/ minval(vdpth) , maxval(vdpth) /)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_dpth,'valid_min',NF_REAL,1,rrange(1)), &
            &        croutnm, cfil, cvar)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_dpth,'valid_max',NF_REAL,1,rrange(2)), &
            &        croutnm, cfil, cvar)
       !!
       !!      For time
       cu = 'unknown'
       IF ( present(cun_t) ) cu = cun_t
       ils = len_trim(trim(cu))
       CALL disp_err(NF_PUT_ATT_TEXT(id_fil, id_time, 'units', ils, trim(cu)), &
            &        croutnm, cfil, cvar)
       rrange = (/ minval(vtime), maxval(vtime) /)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_time,'valid_min',NF_REAL,1,rrange(1)), &
            &        croutnm, cfil, cvar)
       CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_time,'valid_max',NF_REAL,1,rrange(2)), &
            &          croutnm, cfil, cvar)
       !!
       !!
       !! Global attributes
       ils = LEN_TRIM(cabout)
       CALL disp_err(NF_PUT_ATT_TEXT(id_fil, NF_GLOBAL, 'About', ils, trim(cabout)), &
            &        croutnm, cfil, cvar)
       !!
       !!
       !!           END OF DEFINITION
       !!           -----------------
       CALL disp_err(NF_ENDDEF(id_fil), croutnm, cfil, cvar)
       !!
       !!
       !!          WRITE COORDINATES
       !!          -----------------
       !!       Write longitude variable :
       CALL disp_err(NF_PUT_VAR_REAL(id_fil, id_lon, rlon), croutnm, cfil, cvar)
       !!
       !!       Write latitude variable :
       CALL disp_err(NF_PUT_VAR_REAL(id_fil, id_lat, rlat), croutnm, cfil, cvar)
       !!
       !!       Write depth variable :
       CALL disp_err(NF_PUT_VAR_REAL(id_fil, id_dpth, vdpth), croutnm, cfil, cvar)
       !!
       !!       Write time variable :
       CALL disp_err(NF_PUT_VARA_REAL(id_fil, id_time, 1, lt, vtime), &
            &        croutnm, cfil, cvar)
       !! 
      ENDIF   !INITIALIZING FILE


       !!
       !!  INITIALIZING VARIABLES
       !!
       CALL disp_err(NF_REDEF(id_fil), croutnm, cfil, cvar)
       !!
       !!
       !!  DEFINE NEW VARIABLE
       CALL disp_err(NF_DEF_VAR(id_fil, trim(cvar), NF_REAL, 4, id_dims4, id_var), &
              &        croutnm, cfil, cvar)
       !!
       !!  VARIABLE ATTRIBUTES
       !!
       !! Long name
       ils=LEN_TRIM(trim(cln))
       CALL disp_err(NF_PUT_ATT_TEXT(id_fil, id_var, 'long_name', ils, trim(cln)), &
            &       croutnm, cfil, cvar)
       !!
       !! Units
       ils = LEN_TRIM(cunit) 
       CALL disp_err(NF_PUT_ATT_TEXT(id_fil, id_var, 'units', ils, trim(cunit) ),  &
            &       croutnm, cfil, cvar)     
       !!
       !! Missing value
       IF ( vflag /= 0.) &
            & CALL disp_err(NF_PUT_ATT_REAL(id_fil, id_var, 'missing_value', &
            &              NF_REAL, 1, vflag), croutnm, cfil, cvar)
       !!
       !! Calculating min without the mask vflag :
       !IF ( vflag /= 0.) THEN
       !   rmin = 1.E4
       !   DO ji=1, lx
       !      DO jj=1, ly
       !         IF ((x3d(ji,jj,1) <= rmin).and.(x3d(ji,jj,1) > vflag)) &
       !              & rmin = x3d(ji,jj,1)
       !      END DO
       !   END DO
       !ELSE
       !   rmin = minval(x3d)
       !END IF
       !!
       !rrange = (/ rmin , maxval(x3d) /)
       !!
       !CALL disp_err(NF_PUT_ATT_REAL(id_fil,id_var,'valid_range',NF_REAL,2,rrange), &
       !     &        croutnm, cfil, cvar)
       !!
       !!           END OF DEFINITION
       !!           -----------------
       CALL disp_err(NF_ENDDEF(id_fil), croutnm, cfil, cvar)
       
    ELSE    
       !!
       !!   GET IDENTIFIERS NECESSARY FOR ADDING DATA
       !!
       !!       Get Variable Identifier
       CALL disp_err(NF_INQ_VARID(id_fil, trim(cvar), id_var), &
            &        croutnm, cfil, cvar)
       

    END IF  !lct == 1
    !!          
    !!
    !!
    !!          WRITE VARIABLE
    istart4 = (/ 1,  1, 1, lct /)
    icount4 = (/ lx, ly, lz, 1 /)
    CALL disp_err(NF_PUT_VARA_REAL(id_fil, id_var, istart4, icount4, x3d), &
         &        croutnm, cfil, cvar)
    !!
    !!
    IF ( lct == lt ) CALL disp_err(NF_CLOSE(id_fil), croutnm, cfil, cvar)
    !!
    !!
  END SUBROUTINE P3D_T_IRR
  !!
  !!
  !!
  !!
  !!
  !!
  !!
  !!
  SUBROUTINE CHECK_4_MISS(cfil, cvar, lmv, rmissval, cmiss)
    !!
    !! o This routine looks for the presence of a missing value attribute 
    !!   of variable cvar into file cfil
    !!
    !! INPUT :
    !! -------
    !!         * cvar    = variable                                [character]
    !!         * cfil    = treated file                              [character]
    !! 
    !! OUTPUT :
    !! --------
    !!         * imiss    = 0 -> no missing value, 1 -> missing value found   [integer]
    !!         * rmissval = value of missing value                          [real]
    !!         * [cmiss]  = name of the missing value arg. |OPTIONAL|  [character]
    !!
    !! Author : L. BRODEAU, december 2008
    !!
    !!----------------------------------------------------------------------------
    !!
    CHARACTER(len=*), INTENT(in)  :: cfil, cvar
    !! 
    LOGICAL,            INTENT(out) :: lmv
    REAL(wpez),         INTENT(out) :: rmissval
    !!
    CHARACTER(len=20) , OPTIONAL, INTENT(in)  :: cmiss
    !!
    INTEGER :: lx, ly, kz, kt, ierr
    !!
    croutnm = 'CHECK_4_MISS'
    !!
    !!
    !! Opening file :
    CALL disp_err(NF_OPEN(cfil, NF_NOWRITE, nfid), croutnm, cfil, cvar)
    !!
    !! Chosing variable :
    CALL disp_err(NF_INQ_VARID(nfid, cvar, nvid), croutnm, cfil, cvar)
    !!
    !!
    IF ( present(cmiss) ) THEN
       ierr = NF_GET_ATT_REAL(nfid, nvid, cmiss, rmissval)
    ELSE
       !! Default name for a missing value is "missing_value" :
       ierr = NF_GET_ATT_REAL(nfid, nvid, 'missing_value', rmissval)
    END IF
    !!
    IF ( ierr == -43 ) THEN
       lmv = .FALSE.
       !!
    ELSE
       !!
       IF (ierr ==  NF_NOERR) THEN
          lmv = .TRUE.
       ELSE
          PRINT *, 'ERROR getting missing_value attribute: ierr=',ierr
          PRINT *, 'Subroutine CHECK_4_MISS of eZcdf'
          STOP
       END IF
    END IF
    !!
    CALL disp_err(NF_CLOSE(nfid), croutnm, cfil, cvar)
    !!
  END SUBROUTINE CHECK_4_MISS
  !!
  !!
  SUBROUTINE GET_VAR_INFO(cfil, cvar, cunit, clnm)
    !!
    !! o This routine returns the unit and longname of variable if they exist!
    !!
    !! INPUT :
    !! -------
    !!         * cvar    = variable                                [character]
    !!         * cfil    = treated file                            [character]
    !! 
    !! OUTPUT :
    !! --------
    !!         * cunit = unit of cvar                              [character]
    !!         * clnm  = name of the missing value arg.            [character]
    !!
    !! Author : L. BRODEAU, 2008
    !!
    !!----------------------------------------------------------------------------
    !!
    CHARACTER(len=*), INTENT(in)  :: cfil, cvar
    CHARACTER(len=20) , INTENT(out) :: cunit
    CHARACTER(len=200), INTENT(out) :: clnm
    !!
    INTEGER :: lx, ly, kz, kt, ierr
    !!
    croutnm = 'GET_VAR_INFO'
    !!
    !!
    !! Opening file :
    CALL disp_err(NF_OPEN(cfil, NF_NOWRITE, nfid), croutnm, cfil, cvar)
    !!
    !! Chosing variable :
    CALL disp_err(NF_INQ_VARID(nfid, cvar, nvid), croutnm, cfil, cvar)
    !!
    !!
    ierr = NF_GET_ATT_TEXT(nfid, nvid, 'units', cunit)
    IF (ierr /= 0) cunit = 'UNKNOWN'
    !!
    ierr = NF_GET_ATT_TEXT(nfid, nvid, 'long_name', clnm)
    IF (ierr /= 0) clnm = 'UNKNOWN'
    !!
    CALL disp_err(NF_CLOSE(nfid), croutnm, cfil, cvar)
    !!
  END SUBROUTINE GET_VAR_INFO
  !!
  !!
  !!
  SUBROUTINE PRTMASK(lx, ly, xmsk, cfil, cvar,   xlon, xlat, cvarlon, cvarlat)
    !!
    !!-----------------------------------------------------------------------------
    !! This routine prints an integer array of a 2D mask into a netcdf file 
    !! ( healthy minds usually agree for earth == 0, sea == 1 )
    !!
    !! INPUT :
    !! -------
    !!        lx    = x dimension of mask                      [integer]
    !!        ly    = y dimension of mask                      [integer]
    !!        xmsk  = 2D array (lx,ly) contening mask           [integer]
    !!        cfil  = name of the output file                   [character]
    !!        cvar  = name of the mask variable                 [character]
    !! OPTIONAL :
    !! ----------
    !!        xlon    = longitude array
    !!        xlat    = latitude array
    !!        cvarlon = longitude name
    !!        cvarlat = latitude name
    !!
    !!------------------------------------------------------------------------------
    !!
    INTEGER,                      INTENT(in) :: lx, ly
    !!INTEGER,    DIMENSION(lx,ly), INTENT(in) :: xmsk
    REAL(wpez),    DIMENSION(lx,ly), INTENT(in) :: xmsk
    CHARACTER(len=*),             INTENT(in) :: cfil, cvar
    !!
    REAL(wpez), DIMENSION(lx,ly), INTENT(in), OPTIONAL :: xlon, xlat
    CHARACTER(len=*),             INTENT(in), OPTIONAL :: cvarlon, cvarlat

    !!
    LOGICAL :: lzcoord
    !!
    !!
    croutnm = 'PRTMASK'
    !!
    lzcoord = .FALSE.
    !!
    IF ( present(xlon).AND.present(xlat) ) THEN
       IF ( present(cvarlon).AND.present(cvarlat) ) THEN
          lzcoord = .TRUE.
       ELSE
          PRINT *, ''; PRINT *, 'ERROR: If you specify xlon and xlat in PRTMASK (io_ezcdf.f90),'
          PRINT *, 'then you must also specify cvarlon and cvarlat!'; PRINT *, ''
          STOP
       END IF
    END IF
    !!
    !!
    !!           CREATE NETCDF OUTPUT FILE :
    !!           ---------------------------    
    CALL disp_err(NF_CREATE(CFIL,NF_CLOBBER,nfid), croutnm, cfil, cvar)
    !!
    !!
    !!             DIMMENSIONS 
    !!             -----------
    !!       Create longitude dimension :
    CALL disp_err(NF_DEF_DIM(nfid, 'x', lx, id_x), croutnm, cfil, cvar)
    !!
    !!       Create latitude dimension :
    CALL disp_err(NF_DEF_DIM(nfid, 'y', ly, id_y), croutnm, cfil, cvar)
    !!
    !!
    id_dims2 = (/ id_x, id_y /)
    !!
    !!
    !!           VARIABLE TO PLOT
    !!           ----------------
    IF ( lzcoord ) THEN
       !! Longitude
       CALL disp_err(NF_DEF_VAR(nfid, trim(cvarlon), NF_REAL, 2, id_dims2, id_lon), &
            &         croutnm, cfil, cvar)
       !!
       !! Latitude
       CALL disp_err(NF_DEF_VAR(nfid, trim(cvarlat), NF_REAL, 2, id_dims2, id_lat), &
            &         croutnm, cfil, cvar)
    END IF
    !!
    !! Mask
    !CALL disp_err(NF_DEF_VAR(nfid, cvar, NF_INT, 2, id_dims2, nvid), &
    !     &        croutnm, cfil, cvar)
    CALL disp_err(NF_DEF_VAR(nfid, trim(cvar), NF_REAL, 2, id_dims2, nvid), &
         &        croutnm, cfil, cvar)
    !!
    !!
    !!           END OF DEFINITION
    !!           -----------------
    CALL disp_err(NF_ENDDEF(nfid), croutnm, cfil, cvar)
    !!
    !!
    !!          WRITE COORDINATES
    !!          -----------------
    IF ( lzcoord ) THEN
       CALL disp_err(NF_PUT_VAR_REAL(nfid, id_lon, xlon), croutnm, cfil, cvar)
       CALL disp_err(NF_PUT_VAR_REAL(nfid, id_lat, xlat), croutnm, cfil, cvar)
    END IF
    !!
    !!
    !!          WRITE VARIABLE
    !!          --------------
    !CALL disp_err(NF_PUT_VAR_INT(nfid, nvid, xmsk), croutnm, cfil, cvar)
    CALL disp_err(NF_PUT_VAR_REAL(nfid, nvid, xmsk), croutnm, cfil, cvar)
    !NF_PUT_VARA_REAL lolo
    !!
    CALL disp_err(NF_CLOSE(nfid), croutnm, cfil, cvar)
    !!
    !!
  END SUBROUTINE PRTMASK
  !!
  !!
  !!
  SUBROUTINE GET_SF_AO(cfil, cvar, rsf, rao)
    !!
    !!-----------------------------------------------------------------------
    !! This routine extracts the 'scale_factor' and 'add_offset' of a given
    !! variable from a netcdf file
    !!
    !! INPUT :
    !! -------
    !!          * cfil      : name of the input file              (character)
    !!          * cvar      : name of the variable                (character)
    !!
    !! OUTPUT :
    !! --------
    !!          * rsf       : scale factor                        (real)
    !!          * rao       : add offset                          (real)
    !!
    !!
    !!  Author :            Laurent Brodeau, 2004
    !!  --------
    !!------------------------------------------------------------------------
    !!  
    CHARACTER(len=*), INTENT(in)  :: cfil, cvar
    REAL(wpez),         INTENT(out) :: rsf, rao
    !!
    !! local :
    CHARACTER(len=200) :: cname, catt_nm
    INTEGER :: ja, itype, idims, numb_att, jsf, jao
    LOGICAL :: lsf, lao
    !!
    INTEGER, DIMENSION(30) :: vdims
    !!
    !!
    croutnm = 'GET_SF_AO'
    !!
    lsf = .FALSE. ;   lao = .FALSE.
    rsf = 1.      ;   rao = 0.
    !!
    !!
    CALL disp_err(NF_OPEN(cfil, NF_NOWRITE, nfid), croutnm, cfil, cvar)
    !!
    CALL disp_err(NF_INQ_VARID(nfid, cvar, nvid), croutnm, cfil, cvar)
    !!
    CALL disp_err(NF_INQ_VAR(nfid, nvid, cname, itype, idims, vdims, numb_att), &
         &        croutnm, cfil, cvar)
    !!
    !!
    !!
    DO ja = 1, numb_att 
       !!
       CALL disp_err(NF_INQ_ATTNAME(nfid, nvid, ja, catt_nm), croutnm, cfil, cvar)  
       !!
       IF ( trim(catt_nm) == 'scale_factor' ) THEN
          lsf = .TRUE.
          jsf = ja
       END IF
       !!
       IF ( trim(catt_nm) == 'add_offset' ) THEN
          lao = .TRUE.
          jao = ja
       END IF
       !!
    END DO
    !!
    !!
    IF ( lsf .and. lao) THEN
       !!
       !! scale-factor :
       CALL disp_err(NF_GET_ATT_REAL(nfid, nvid, 'scale_factor', rsf), &
            &        croutnm, cfil, cvar)
       !!
       !! Add-offset :
       CALL disp_err(NF_GET_ATT_REAL(nfid, nvid, 'add_offset',   rao), &
            &        croutnm, cfil, cvar)
       !!
       !! Closing file :
       CALL disp_err(NF_CLOSE(nfid), croutnm, cfil, cvar)
       !!
       !!
    END IF
    !!
  END SUBROUTINE GET_SF_AO
  !!
  !!
  !!
  !!
  !!===============================================================================
  !!
  !! E N D
END MODULE io_ezcdf
!!
!!
