      subroutine vfthil(n, no, x, stc)
c routine written by Richard James for 3-D Poisson solver
      integer n, no
      logical stc
      real x( n * no )
c
c     purpose
c
c        to perform a hilbert transformation on a 2-dimensional
c        mesh held in a 1-dimensional array.  the transform direction
c        is orthogonal to the store direction.
c
c     parameters
c
c     n    =  transform length
c
c     no   =  number of transforms to be evaluated in parallel
c
c     x    -  a real array holding the data
c
c     stc  -  (type logical) - requests that the hilbert transform
c             be used to change a sine transform into a cosine
c             transform.
c
c     note:  the routine aborts if the transform length exceeds 257.
c            this limit may be increased by changing the size of the
c            local integer array  w.  vectors are assumed to be
c            contiguous in memory, and so orthogonal to the transform
c            direction.
c
      include 'rjinclude.h'
c
c local variables
      integer i, iw, iw2, i0, i2, i3, i4, j, jj, k, nb, nm1
      real recip, u, v, xx1, xx2
c
c      logical stc
c      real x(n*no)
      integer st0, st1, st2, st3, st4, st5, set, step
      integer xx, y1, y2, z1, z2, sum, diff, lw(192)
c
c     set initial values and check transform length
c
      nm1 = n - 1
      recip = 1.0/float(nm1)
      nb = nm1/2
      step = no
      if(nb.gt.128) then
        write( s2, 200 )n
 200    format( 'Increase pointer array size in vfthil(xd)' /
     1          ' transform length =', i10/
     2          ' we need 3*(transform length - 1)/4 pointers' )
        call crash( 'vfthil', 'xd 1')
      end if
c
c     set up pointers for work areas in /lev3/.
c
      xx = 0
      y1 = xx + step
      y2 = y1 + step
      z1 = y2 + step
      z2 = z1 + step
      sum = z2 + step
      diff = sum + step
c
c     check space in /lev3/ is adequate.
c
      i = diff + step
      if(i.gt.lnlev3) then
        write( s2, 201 )i, lnlev3
 201    format( 'Insufficient space in /lev3/ in vfthil(xd)' /
     1          ' space required is', i10, ' words,',
     2          ' space available is', i10, ' words' )
        call crash( 'vfthil', 'xd 2')
        stop
      end if
c
c     set up pointers for work areas in stack in blank common.
c
      do 1 i = 0, 3*nb/2 - 1
 1    lw(i + 1) = istack + i*step
c
c     check for stack overflow.
c
      j = step*max0(3*nb/2, 5) + istack
      if(j.gt.maxstk) then
        maxstk = j
        mxstid = 'vfthil'
      end if
      if(j.gt.lstack) then
        write( s2, 202 )istack, j, lstack
 202    format( 'Stack overflow in vfthil(xd)' /
     1          ' stack pointer =', i10, ' last word required =',
     2          i10, ' limit =', i10 )
        call crash( 'vfthil', 'xd 3')
        stop
      end if
c
c     select prologue according to transform required
c
      if(stc) go to 3
c
c     prepare cos-to-sine transformation
c
      st2 = nb*step
      do 50 i = 1, step
50    lev3(xx + i) = x(st2 + i)
      do 51 i = 1, step
51    lev3(y1 + i) = 0.0
c
c     compute end results
c
      st1 = - step
      do 2 j = 1, nb
      st1 = st1 + step
      st2 = st2 + step
      do 52 i = 1, step
52    lev3(xx + i) = lev3(xx + i) + x(st1 + i)
      do 53 i = 1, step
53    lev3(y1 + i) = lev3(y1 + i) + x(st2 + i)
 2    continue
      do 54 i = 1, step
54    vect(i) = recip*(lev3(xx + i) + lev3(y1 + i))
      do 55 i = 1, step
55    lev3(sum + i) = vect(i)
      do 56 i = 1, step
56    vect(i) = recip*(lev3(xx + i) - lev3(y1 + i))
      do 57 i = 1, step
57    lev3(diff + i) = vect(i)
      go to 5
c
c     prepare sine-to-cos transformation
c
c     preserve sum and difference of end values
c
 3    st0 = nm1*step
      do 58 i = 1, step
58    lev3(sum + i) = 2.0*(x(i) + x(st0 + i))
      do 59 i = 1, step
59    lev3(diff + i) = 2.0*(x(i) - x(st0 + i))
c
c     initial pointer for twiddle factors
c
      iw = nb
c
c     cycle over data for adjustment
c
      do 4 j = 9, n, 4
      st1 = st0 - step
      st2 = st1 - step
      st3 = st2 - step
      u = fact(iw - 1)
      v = fact(iw)
      iw = iw - 2
c
c     adjust first five sets of values in current set
c
      i0 = st0
      st0 = st3 - step
      do 60 i = 1, step
60    vect(i) = u*x(st1 + i)
      do 61 i = 1, step
61    x(i0 + i) = vect(i)
      do 62 i = 1, step
62    vect(i) = u*x(st2 + i)
      do 63 i = 1, step
63    x(st1 + i) = vect(i)
      do 64 i = 1, step
64    vect(i) = v*x(st3 + i)
      do 65 i = 1, step
65    x(st2 + i) = vect(i)
      do 66 i = 1, step
66    vect(i) = v*x(st0 + i)
      do 67 i = 1, step
67    x(st3 + i) = vect(i)
 4    continue
c
c     adjust first five values in first set
c
      st1 = st0 - step
      st2 = st1 - step
      st3 = st2 - step
      u = fact(2)
      do 68 i = 1, step
68    vect(i) = u*x(st1 + i)
      do 69 i = 1, step
69    x(st0 + i) = vect(i)
      do 70 i = 1, step
70    vect(i) = u*x(st2 + i)
      do 71 i = 1, step
71    x(st1 + i) = vect(i)
      do 72 i = 1, step
72    vect(i) = x(st3 + i)
      do 73 i = 1, step
73    x(st2 + i) = vect(i)
      do 74 i = 1, step
74    x(st3 + i) = 0.0
      do 75 i = 1, step
75    x(i)  = 0.0
c
c     evaluation of s-alphas
c
c     calculate t-alphas from odd data
c
 5    st1 = step*(nb + 1)
      do 6 jj = 1, nb
      do 76 i = 1, step
76    wstack(lw(jj) + i) = x(st1 + i)
 6    st1 = st1 + step
c
c     initial number of sets for folding
c
      set = 1
c
c     advance to next fold
c
 7    nb = nb/4
      st1 = 1
      if((16*set + 1).gt.n) go to 10
      iw = nb + 1
      iw2 = iw + nb
c
c     cycle over sets for folding
c
      do 9 j = 1, nb
c
c     step pointers
c
      st2 = st1 + set
      st3 = st2 + set
      st4 = st3 + set
c
c     find twiddle factors
c
      u = fact(iw)
      xx1 = fact(iw2)
      xx2 = fact(iw2 + 1)
      iw = iw + 1
      iw2 = iw2 + 2
c
c     fold current set
c
      do 8 k = 1, set
      do 77 i = 1, step
77    lev3(y1 + i) = xx1*(wstack(lw(st1) + i) - wstack(lw(st2) + i))
      do 78 i = 1, step
78    lev3(z1 + i) = wstack(lw(st1) + i) + wstack(lw(st2) + i)
      do 79 i = 1, step
79    lev3(y2 + i) = xx2*(wstack(lw(st3) + i) - wstack(lw(st4) + i))
      do 80 i = 1, step
80    lev3(z2 + i) = wstack(lw(st3) + i) + wstack(lw(st4) + i)
      do 81 i = 1, step
81    wstack(lw(st1) + i) = u*(lev3(y1 + i) - lev3(y2 + i))
      do 82 i = 1, step
82    wstack(lw(st3) + i) = lev3(y1 + i) + lev3(y2 + i)
      do 83 i = 1, step
83    wstack(lw(st2) + i) = u*(lev3(z1 + i) - lev3(z2 + i))
      do 84 i = 1, step
84    wstack(lw(st4) + i) = lev3(z1 + i) + lev3(z2 + i)
      st1 = st1 + 1
      st2 = st2 + 1
      st3 = st3 + 1
 8    st4 = st4 + 1
 9    st1 = st4
c
c     increase set length and advance to next fold
c
      set = 4*set
      go to 7
c
c     check for single fold
c
10    if((8*set + 1).gt.n) go to 13
c
c     perform single fold
c
      st2 = set + 1
      u = fact(3)
      do 12 j = 1, 2
      do 11 k = 1, set
      do 85 i = 1, step
85    lev3(y1 + i) = u*(wstack(lw(st1) + i) - wstack(lw(st2) + i))
      do 86 i = 1, step
86    vect(i) = wstack(lw(st1) + i) + wstack(lw(st2) + i)
      do 87 i = 1, step
87    wstack(lw(st2) + i) = vect(i)
      do 88 i = 1, step
88    wstack(lw(st1) + i) = lev3(y1 + i)
      st1 = st1 + 1
11    st2 = st2 + 1
      st1 = st2
      st2 = st1 + set
12    u = fact(4)
c
c     initialise counts for unfolding
c
13    set = nm1/4
      nb = 1
c
c     generate lowest level t-alphas for unfolding
c
      st1 = 1
      st2 = set + 1
      st3 = set + st2
      u = fact(2)
      do 14 j = 1, set
      do 89 i = 1, step
89    vect(i) = u*(wstack(lw(st2) + i) - wstack(lw(st1) + i))
      do 90 i = 1, step
90    wstack(lw(st3) + i) = vect(i)
      do 91 i = 1, step
91    vect(i) = wstack(lw(st2) + i) + wstack(lw(st1) + i)
      do 92 i = 1, step
92    wstack(lw(st2) + i) = vect(i)
      do 93 i = 1, step
93    vect(i) = wstack(lw(st3) + i)
      do 94 i = 1, step
94    wstack(lw(st1) + i) = - vect(i)
      st1 = st1 + 1
      st2 = st2 + 1
14    st3 = st3 + 1
c
c     initial pointers for unfolding process
c
      st0 = nm1/2 + 1
      set = set/2
c
c     initial pointers for current unfolding stage
c
15    st1 = 1
      st2 = set + 1
      st3 = set + st2
      st4 = set + st3
      st5 = st0 - set
c
c     unfold first quartet
c
      do 16 j = 1, set
      do 95 i = 1, step
95    lev3(y1 + i) = wstack(lw(st1) + i) + wstack(lw(st3) + i)
      do 96 i = 1, step
96    lev3(y2 + i) = wstack(lw(st2) + i) + wstack(lw(st4) + i)
      do 97 i = 1, step
97    vect(i) = wstack(lw(st1) + i) - wstack(lw(st3) + i)
      do 98 i = 1, step
98    wstack(lw(st5) + i) = vect(i)
      do 99 i = 1, step
99    wstack(lw(st1) + i) = lev3(y1 + i) + lev3(y2 + i)
      do 100 i = 1, step
100   wstack(lw(st2) + i) = lev3(y1 + i) - lev3(y2 + i)
      st1 = st1 + 1
      st2 = st2 + 1
      st3 = st3 + 1
      st4 = st4 + 1
16    st5 = st5 + 1
c
c     adjust zero and pack odd section at end
c
      st0 = st0 - set - (set + 1)/2
      st2 = st0
      do 17 j = 1, set, 2
      do 101 i = 1, step
101   vect(i) = wstack(lw(st1) + i)
      do 102 i = 1, step
102   wstack(lw(st2) + i) = vect(i)
      st1 = st1 + 1
17    st2 = st2 + 1
c
c     next set requires single fold only
c
      st1 = st2 + set
      st2 = st1 + set
      u = fact(2)
      do 18 j = 1, set
      do 103 i = 1, step
103   lev3(y1 + i) = u*wstack(lw(st2) + i)
      do 104 i = 1, step
104   vect(i) = wstack(lw(st1) + i) - lev3(y1 + i)
      do 105 i = 1, step
105   wstack(lw(st2) + i) = vect(i)
      do 106 i = 1, step
106   wstack(lw(st1) + i) = wstack(lw(st1) + i) + lev3(y1 + i)
      st1 = st1 + 1
18    st2 = st2 + 1
c
c     jump if initial unfolding
c
      if(nb.eq.1) go to 21
c
c     initial pointers
c
      iw = 2
      iw2 = 3
      st1 = st2
c
c     cycle over sets for unfolding
c
      do 20 jj = 2, nb
c
c     find twiddle factors
c
      u = fact(iw)
      xx1 = fact(iw2)
      xx2 = fact(iw2 + 1)
      iw = iw + 1
      iw2 = iw2 + 2
c
c     initial pointers
c
      st2 = st1 + set
      st3 = st2 + set
      st4 = st3 + set
c
c     fold quartet
c
      do 19 j = 1, set
      do 107 i = 1, step
107   lev3(y1 + i) = wstack(lw(st1) + i) + u*wstack(lw(st3) + i)
      do 108 i = 1, step
108   lev3(z1 + i) = wstack(lw(st1) + i) - u*wstack(lw(st3) + i)
      do 109 i = 1, step
109   lev3(y2 + i) = wstack(lw(st2) + i) + u*wstack(lw(st4) + i)
      do 110 i = 1, step
110   lev3(z2 + i) = wstack(lw(st2) + i) - u*wstack(lw(st4) + i)
      do 111 i = 1, step
111   wstack(lw(st1) + i) = lev3(y1 + i) + xx1*lev3(y2 + i)
      do 112 i = 1, step
112   wstack(lw(st2) + i) = lev3(y1 + i) - xx1*lev3(y2 + i)
      do 113 i = 1, step
113   wstack(lw(st3) + i) = lev3(z1 + i) + xx2*lev3(z2 + i)
      do 114 i = 1, step
114   wstack(lw(st4) + i) = lev3(z1 + i) - xx2*lev3(z2 + i)
      st1 = st1 + 1
      st2 = st2 + 1
      st3 = st3 + 1
19    st4 = st4 + 1
c
c     pointer for next set
c
20    st1 = st4
c
c     test for completion
c
21    if(set.eq.1) go to 24
      if(set.eq.2) go to 22
      set = set/4
      nb = 4*nb
      go to 15
c
c     single fold
c
22    st0 = st0 - 2
      do 115 i = 1, step
115   vect(i) = wstack(lw(1) + i) + wstack(lw(2) + i)
      do 116 i = 1, step
116   wstack(lw(st0) + i) = vect(i)
      do 117 i = 1, step
117   vect(i) = wstack(lw(1) + i) - wstack(lw(2) + i)
      do 118 i = 1, step
118   wstack(lw(st0 + 1) + i) = vect(i)
      iw = 2
      st1 = st0 + 3
      do 23 j = 9, n, 4
      u = fact(iw)
      do 119 i = 1, step
119   lev3(y1 + i) = u*wstack(lw(st1 + 1) + i)
      do 120 i = 1, step
120   vect(i) = wstack(lw(st1) + i) - lev3(y1 + i)
      do 121 i = 1, step
121   wstack(lw(st1 + 1) + i) = vect(i)
      do 122 i = 1, step
122   wstack(lw(st1) + i) = wstack(lw(st1) + i) + lev3(y1 + i)
      iw = iw + 1
23    st1 = st1 + 2
      go to 25
24    st0 = st0 - 1
      do 150 i = 1, step
150   vect(i) = wstack(lw(1) + i)
      do 153 i = 1, step
153   wstack(lw(st0) + i) = vect(i)
c
c     return t-alphas to result area and pick up data for u-alphas.
c     this is stacked in the working area with blanks after the first
c     three elements, to allow for shifting down.
c
25    nb = nm1/2 + 1
      st1 = step*nb - step
      st3 = 3*nm1/4
      do 26 j = 1, nb
      do 123 i = 1, step
123   lev3(y1 + i) = x(st1 + i)
      do 124 i = 1, step
124   x(st1 + i) = wstack(lw(st3) + i)
      do 125 i = 1, step
125   wstack(lw(st3) + i) = lev3(y1 + i)
      st1 = st1 - step
26    st3 = st3 - 1
      do 126 i = 1, step
126   vect(i) = wstack(lw(st3 + 1) + i)
      do 127 i = 1, step
127   wstack(lw(1) + i) = vect(i)
      do 128 i = 1, step
128   vect(i) = wstack(lw(st3 + 2) + i)
      do 129 i = 1, step
129   wstack(lw(2) + i) = vect(i)
      do 130 i = 1, step
130   vect(i) = wstack(lw(st3 + 3) + i)
      do 131 i = 1, step
131   wstack(lw(4) + i) = vect(i)
c
c     initialise folding of even data
c
      set = 1
      nb = nm1/8
      st0 = 6*nb + 1
      if(nb.lt.2) go to 33
c
c     fold to next level
c
27    st1 = st0
      st2 = st0 + set
      st3 = st2 + set
      st4 = st3 + set
      iw = nb + 1
      iw2 = iw + iw
c
c     cycle over sets for folding
c
      do 28 jj = 2, nb
c
c     step pointers
c
      i3 = 3*set
      st1 = st1 - i3
      st2 = st2 - i3
      st3 = st3 - i3
      st4 = st4 - i3
      iw = iw - 1
      iw2 = iw2 - 2
c
c     pick up twiddle factors
c
      u = fact(iw)
      xx1 = fact(iw2 - 1)
      xx2 = fact(iw2)
c
c     fold set
c
      do 28 j = 1, set
      st1 = st1 - 1
      st2 = st2 - 1
      st3 = st3 - 1
      st4 = st4 - 1
      do 132 i = 1, step
132   lev3(y1 + i) = xx1*(wstack(lw(st1) + i) - wstack(lw(st2) + i))
      do 133 i = 1, step
133   lev3(z1 + i) = wstack(lw(st1) + i) + wstack(lw(st2) + i)
      do 134 i = 1, step
134   lev3(y2 + i) = xx2*(wstack(lw(st3) + i) - wstack(lw(st4) + i))
      do 135 i = 1, step
135   lev3(z2 + i) = wstack(lw(st3) + i) + wstack(lw(st4) + i)
      do 136 i = 1, step
136   wstack(lw(st1) + i) = u*(lev3(y1 + i) - lev3(y2 + i))
      do 137 i = 1, step
137   wstack(lw(st3) + i) = lev3(y1 + i) + lev3(y2 + i)
      do 138 i = 1, step
138   wstack(lw(st2) + i) = u*(lev3(z1 + i) - lev3(z2 + i))
      do 139 i = 1, step
139   wstack(lw(st4) + i) = lev3(z1 + i) + lev3(z2 + i)
28    continue
c
c     fold top set and move back
c
      i2 = set + set
      i4 = i2 + i2
      st1 = st1 - i2
      st2 = st1 + set
      st5 = i4 + 1
      st3 = st5 + i2
      st4 = st3 + set
      u = fact(2)
      do 29 j = 1, set
      do 140 i = 1, step
      vect(i) = wstack(lw(st1) + i) + wstack(lw(st2) + i)
140   lev3(y1 + i) = wstack(lw(st1) + i) - wstack(lw(st2) + i)
      do 141 i = 1, step
141   wstack(lw(st4) + i) = vect(i)
      do 142 i = 1, step
142   wstack(lw(st3) + i) = u*lev3(y1 + i)
      do 143 i = 1, step
143   wstack(lw(st5) + i) = 0.0
      do 144 i = 1, step
144   wstack(lw(st5 + 1) + i) = 0.0
      st1 = st1 + 1
      st2 = st2 + 1
      st3 = st3 + 1
      st4 = st4 + 1
29    st5 = st5 + 2
c
c     move next set
c
      if(nb.eq.2) go to 31
      st1 = st2
      st2 = st4 + i4
      do 30 j = 1, i4
      do 145 i = 1, step
145   vect(i) = wstack(lw(st1) + i)
      do 146 i = 1, step
146   wstack(lw(st2) + i) = vect(i)
      do 147 i = 1, step
147   wstack(lw(st4) + i) = 0.0
      st1 = st1 + 1
      st2 = st2 + 1
30    st4 = st4 + 1
c
c     fold first quartet
c
31    st1 = 1
      st2 = st1 + set
      st3 = st2 + set
      st4 = st3 + set
      do 32 j = 1, set
      do 148 i = 1, step
      vect(i) = wstack(lw(st1) + i) - wstack(lw(st2) + i)
148   lev3(z1 + i) = wstack(lw(st1) + i) + wstack(lw(st2) + i)
      do 149 i = 1, step
149   wstack(lw(st1) + i) = vect(i)
      do 151 i = 1, step
151   wstack(lw(st3) + i) = vect(i)
      do 152 i = 1, step
152   vect(i) = wstack(lw(st4) + i)
      do 154 i = 1, step
154   wstack(lw(st2) + i) = lev3(z1 + i) - vect(i)
      do 155 i = 1, step
155   wstack(lw(st4) + i) = lev3(z1 + i) + vect(i)
      st1 = st1 + 1
      st2 = st2 + 1
      st3 = st3 + 1
32    st4 = st4 + 1
      set = i4
c
c     test for folding complete
c
      if(nb.eq.2) go to 37
      if(nb.eq.4) go to 33
c
c     advance to next fold
c
      nb = nb/4
      go to 27
c
c     single fold
c
33    st1 = st0 - set
      st2 = st0
      u = fact(2)
c
c     fold set
c
      do 34 j = 1, set
      st1 = st1 - 1
      st2 = st2 - 1
      do 156 i = 1, step
      vect(i) = wstack(lw(st1) + i) + wstack(lw(st2) + i)
156   lev3(y1 + i) = u*(wstack(lw(st1) + i) - wstack(lw(st2) + i))
      do 157 i = 1, step
157   wstack(lw(st2) + i) = vect(i)
      do 158 i = 1, step
158   wstack(lw(st1) + i) = lev3(y1 + i)
34    continue
c
c     clear zero multiplier locations
c
      i2 = set + set
      st2 = i2 + 1
      do 35 j = 1, set
      do 159 i = 1, step
159   wstack(lw(st2) + i) = 0.0
35    st2 = st2 + 1
c
c     fold first set
c
      st1 = 1
      st2 = st1 + set
      do 36 j = 1, set
      do 160 i = 1, step
      vect(i) = wstack(lw(st1) + i) + wstack(lw(st2) + i)
160   lev3(y1 + i) = wstack(lw(st1) + i) - wstack(lw(st2) + i)
      do 161 i = 1, step
161   wstack(lw(st2) + i) = vect(i)
      do 162 i = 1, step
162   wstack(lw(st1) + i) = lev3(y1 + i)
      st1 = st1 + 1
36    st2 = st2 + 1
      set = i2
c
c     calculate lowest level u-alphas
c
37    st1 = 1
      st2 = set + 1
      st3 = set + st2
      do 38 j = 1, set
      do 163 i = 1, step
      vect(i) = wstack(lw(st2) + i) - wstack(lw(st1) + i)
163   lev3(xx + i) = wstack(lw(st1) + i) + wstack(lw(st2) + i)
      do 165 i = 1, step
165   wstack(lw(st1) + i) = vect(i)
      do 166 i = 1, step
166   vect(i) = wstack(lw(st3) + i) - lev3(xx + i)
      do 167 i = 1, step
167   wstack(lw(st2) + i) = vect(i)
      st1 = st1 + 1
      st2 = st2 + 1
38    st3 = st3 + 1
c
c     initialise unfolding
c
      nb = 1
      iw = 1
c
c     jump if one fold only
c
      if(set.eq.1) go to 41
c
c     initial pointers
c
39    st5 = set/2
      st4 = 1
      st3 = 1 - st5
      st2 = st3 - st5
      st1 = st2 - st5
      iw = nb
      iw2 = nb + nb - 1
c
c     unfold two levels
c
      i3 = set + st5
      do 40 jj = 1, nb
c
c     advance pointers
c
      st1 = st1 + i3
      st2 = st2 + i3
      st3 = st3 + i3
      st4 = st4 + i3
c
c     pick up twiddle factors
c
      iw = iw + 1
      u = fact(iw)
      iw2 = iw2 + 2
      xx1 = fact(iw2)
      xx2 = fact(iw2 + 1)
c
c     unfold quartet of sets
c
      do 40 j = 1, st5
      do 168 i = 1, step
168   lev3(y1 + i) = wstack(lw(st1) + i) + u*wstack(lw(st3) + i)
      do 169 i = 1, step
169   lev3(y2 + i) = wstack(lw(st1) + i) - u*wstack(lw(st3) + i)
      do 170 i = 1, step
170   lev3(z1 + i) = wstack(lw(st2) + i) + u*wstack(lw(st4) + i)
      do 171 i = 1, step
171   lev3(z2 + i) = wstack(lw(st2) + i) - u*wstack(lw(st4) + i)
      do 172 i = 1, step
172   wstack(lw(st1) + i) = lev3(y1 + i) + xx1*lev3(z1 + i)
      do 173 i = 1, step
173   wstack(lw(st2) + i) = lev3(y1 + i) - xx1*lev3(z1 + i)
      do 174 i = 1, step
174   wstack(lw(st3) + i) = lev3(y2 + i) + xx2*lev3(z2 + i)
      do 175 i = 1, step
175   wstack(lw(st4) + i) = lev3(y2 + i) - xx2*lev3(z2 + i)
      st1 = st1 + 1
      st2 = st2 + 1
      st3 = st3 + 1
40    st4 = st4 + 1
c
c     test for completion
c
      if(set.eq.2) go to 43
      if(set.eq.4) go to 41
c
c     advance to next fold
c
      nb = 4*nb
      set = set/4
      go to 39
c
c     single fold
c
41    st1 = - 1
      iw = nm1/4
      nb = iw
      do 42 j = 1, nb
      st1 = st1 + 2
      iw = iw + 1
      u = fact(iw)
      do 176 i = 1, step
      vect(i) = wstack(lw(st1) + i)
176   lev3(y1 + i) = u*wstack(lw(st1 + 1) + i)
      do 177 i = 1, step
177   wstack(lw(st1 + 1) + i) = vect(i) - lev3(y1 + i)
      do 178 i = 1, step
178   wstack(lw(st1) + i) = wstack(lw(st1) + i) + lev3(y1 + i)
42    continue
43    nb = nm1/2
c
c     transfer back to data area
c
      st1 = nb*step + step
      st2 = 1
      do 44 j = 1, nb
      do 179 i = 1, step
179   x(st1 + i) = wstack(lw(st2) + i)
      st1 = st1 + step
44    st2 = st2 + 1
c
c     epilogue
c
      if(stc) go to 46
c
c     finish calculation of sine transform
c
      st1 = step
      st2 = st1 + step
      u = - fact(2)
      st3 = st2 + step
      st0 = st3 + step
      do 180 i = 1, step
180   x(i) = lev3(sum + i)
      do 181 i = 1, step
181   vect(i) = - x(st2 + i)
      do 182 i = 1, step
182   x(st1 + i) = vect(i)
      do 183 i = 1, step
183   vect(i) = u*x(st3 + i)
      do 184 i = 1, step
184   x(st2 + i) = vect(i)
      do 185 i = 1, step
185   vect(i) = u*x(st0 + i)
      do 186 i = 1, step
186   x(st3 + i) = vect(i)
c
c     cycle over rest of data
c
      iw = 3
      do 45 j = 9, n, 4
      st1 = st0 + step
      st2 = st1 + step
      st3 = st2 + step
      i0 = st0
      st0 = st3 + step
      u = - fact(iw + 1)
      v = - fact(iw)
      iw = iw + 2
      do 187 i = 1, step
187   vect(i) = u*x(st1 + i)
      do 188 i = 1, step
188   x(i0 + i) = vect(i)
      do 189 i = 1, step
189   vect(i) = u*x(st2 + i)
      do 190 i = 1, step
190   x(st1 + i) = vect(i)
      do 191 i = 1, step
191   vect(i) = v*x(st3 + i)
      do 192 i = 1, step
192   x(st2 + i) = vect(i)
      do 193 i = 1, step
193   vect(i) = v*x(st0 + i)
      do 194 i = 1, step
194   x(st3 + i) = vect(i)
45    continue
      do 195 i = 1, step
195   x(st0 + i) = lev3(diff + i)
      go to 48
c
c     complete calculation of cosine transform
c
46    st1 = 0
      st2 = nm1*step/2 + step
      do 47 j = 3, n, 2
      do 196 i = 1, step
196   x(st1 + i) = x(st1 + i) + lev3(sum + i)
      do 197 i = 1, step
197   x(st2 + i) = x(st2 + i) + lev3(diff + i)
      st1 = st1 + step
47    st2 = st2 + step
      do 198 i = 1, step
198   x(st1 + i) = x(st1 + i) + lev3(sum + i)
c
c     finish
c
48    continue
      return
      end
