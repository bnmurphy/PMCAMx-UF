      subroutine rndnum(seed,randx)
c
      real randx
c
      integer seed
c
c
c CREATE RANDOM NUMBER
c
c
      seed=2045*seed+1
      seed=seed-(seed/1048576)*1048576
      randx=real(seed+1)/1048577.0
c
c
      return
      end
