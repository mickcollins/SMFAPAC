      subroutine piksrt(n,arr)
      integer arr(n)
      do 12 j=2,n
        a=arr(j)
        do 11 i=j-1,1,-1
         if(arr(i).le.a) go to 10
         arr(i+1)=arr(i)
11      continue
       i=0
10    arr(i+1)=a
12    continue
      return
      end
