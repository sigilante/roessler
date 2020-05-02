/+  roessler, lazytrig
!:
::  3D Projection of Rössler Strange Attractor System
::  
::  By ~lagrev-nocfep (N E Davis)
::
::  Mathematical implementation based on
::    https://en.wikipedia.org/wiki/3D_projection
::
::  Usage:
::    *%/web/output/txt +proj-3d [100 10.000]
::
|=  [ndim=@ud np=@ud]
=/  xmin=@rs  .-0
=/  xmax=@rs  .15
=/  ymin=@rs  .-30
=/  ymax=@rs  .20
=/  dx  (div:rs (sub:rs xmax xmin) (sun:rs ndim))
=/  dy  (div:rs (sub:rs ymax ymin) (sun:rs ndim))
::=<  (gen3d-ppm ndim np)
::=<  (gen3d np)
::=<  (genxyzs:roessler np .0.01)
=<  (gen3d-ppm (ravel-3d-pts (gen-3d-pts ndim np) ndim))
|%
++  gen3d-ppm
  ::  Produce a PPM-styled output.
  |=  [grid=(list (list @rs))]  ^-  (list @t)
  =/  ndim  (lent (snag 0 grid))
  =/  header  ~['P3' (crip "{<ndim>} {<ndim>}") '255']
  =/  ppm-3d-pts  (gen-3d-pts-tape grid)
  (weld header ppm-3d-pts)
++  gen-3d-pts-tape
  ::  Produce a PPM-compatible set of rows.
  |=  [grid=(list (list @rs))]  ^-  (list @t)
  =/  ndim=@ud  (lent (snag 0 grid))
  =/  j=@ud  0
  =/  acc=(list @t)  ~
  ~&  grid
  |-  ^-  (list @t)
    ?:  =(j ndim)  acc
    ~&  "Row: {<j>} of {<ndim>}"
    =/  row  (snag j grid)
    =/  row-text  `@t`(crip (gen-3d-pts-row row))
  $(j +(j), acc `(list @t)`[row-text acc])
++  gen-3d-pts-row
  |=  [row=(list @rs)]  ^-  tape
  =/  ndim=@ud  (lent row)
  =/  i=@ud  0
  =/  acc=tape  ""
  |-  ^-  tape
    ?:  =(i ndim)  acc
    =/  value=@rs  (snag i row)
  $(i +(i), acc (weld "{<value>} {<value>} {<value>} " acc))
++  ravel-3d-pts
  |=  [vec=(list @rs) ndim=@ud]  ^-  (list (list @rs))
  ::%-  flop
  =/  i=@ud  1
  =/  j=@ud  1
  =/  grid=(list (list @rs))  ~
  |-  ^-  (list (list @rs))
    ?:  =(j +(ndim))  grid
    =/  row-start  (mul (sub j 1) ndim)
    =/  row-end    (mul j ndim)
    =/  row  (slag row-start (scag row-end vec))
  $(j +(j), grid `(list (list @rs))`[row grid])
++  gen-3d-pts
  ::  Produce a grid of points for the canvas with rendered points indicated
  |=  [ndim=@ud np=@ud]  ^-  (list @rs)
  =/  pts  (gen3d np)
  ::=/  sorted-pts  (sort-2d-points pts)
  =/  grid  (gengrid ndim)
  =/  index  0
  |-  ^-  (list @rs)
    ?:  =(index np)  grid
    =/  pt=(list @rs)  (snag index pts) ::sorted-pts)
    =/  sx=@rs  i.-:pt
    =/  sy=@rs  i.+<:pt
    ::  Locate the pi-th bin from the floating-point coordinates of the point.
    =/  px=@ud  (bin sx xmin dx np)
    =/  py=@ud  (bin sy ymin dy np)
    =/  pi=@ud  (add px (mul ndim py))
    ::~&  "({<sx>},{<sy>}) -> {<px>} {<py>}, {<pi>}"
    ~&  "Point #{<index>}"
  $(index +(index), grid `(list @rs)`(nick-rs [pi .1] grid))
++  gengrid
  ::  Let's set the grid up as a linear list of point values for now and use
  ::  modulus-style ops to access points.
  |=  [ndim=@ud]  ^-  (list @rs)
  (reap (mul ndim ndim) .0)
++  bin
  ::  Bin the value.
  |=  [p=@rs min=@rs ds=@rs np=@ud]
  =/  index=@ud  0
  |-  ^-  @ud
    ::  The naive algorithm is to count up in range until we find the right one.
    =/  current-bin  (add:rs min (mul:rs ds (sun:rs index)))
    ?:  (lth:rs p current-bin)  index
  $(index +(index))
++  nick-rs
  ::  Replace element of c at index a with item b
  |=  [[a=@ b=@rs] c=(list @rs)]  ^-  (list @rs)
  ?:  =(b 0)  `(list @rs)`c
  (weld (scag a `(list @rs)`c) [b (slag +(a) `(list @rs)`c)])
++  sort-2d-points
  ::  Sort 2D points based on y-value.  Since these are floats, we anticipate
  ::  no collisions and don't sort along the x-axis at all.
  |=  [n=(list (list @rs))]  ^-  (list (list @rs))
  (sort n sort-2d-point-order)
++  sort-2d-point-order
  ::  Compare two 2D points based on their y-axis.
  |=  [a=(list @rs) b=(list @rs)]  ^-  ?
  =/  ay  i.+>-:a
  =/  by  i.+>-:b
  (gth ay by)
++  gen3d
  |=  [np=@ud]  ^-  (list (list @rs))
  =/  pts  (genxyz-pts:roessler np .0.01)
  =/  cx  .114
  =/  cy  .-83
  =/  cz  .100
  =/  theta-x  .0.0
  =/  theta-y  .0.0
  =/  theta-z  .0.0
  =/  ex  .10
  =/  ey  .10
  =/  ez  .10
  =/  index  0
  =/  acc=(list (list @rs))  ~
  |-  ^-  (list (list @rs))
    ?:  =(index np)  acc
    =/  pt  (snag index pts)
    =/  at  i.-:pt
    =/  ax  i.+<:pt
    =/  ay  i.+>-:pt
    =/  az  i.+>+<:pt
    =/  vt  (render-3d-point ax ay az cx cy cz theta-x theta-y theta-z ex ey ez)
  $(index +(index), acc `(list (list @rs))`[vt acc])
++  render-3d-point
  |=  [ax=@rs ay=@rs az=@rs cx=@rs cy=@rs cz=@rs theta-x=@rs theta-y=@rs theta-z=@rs ex=@rs ey=@rs ez=@rs]  ^-  (list @rs)
  =/  sin-x  (sine:lazytrig theta-x)
  =/  sin-y  (sine:lazytrig theta-y)
  =/  sin-z  (sine:lazytrig theta-z)
  =/  cos-x  (cosine:lazytrig theta-x)
  =/  cos-y  (cosine:lazytrig theta-y)
  =/  cos-z  (cosine:lazytrig theta-z)
  =/  xx  (sub:rs ax cx)
  =/  yy  (sub:rs ay cy)
  =/  zz  (sub:rs az cz)
  =/  dx  (sub:rs (mul:rs cos-y (add:rs (mul:rs sin-z yy) (mul:rs cos-z xx))) (mul:rs sin-y zz))
  =/  dy  (add:rs (mul:rs sin-x (add:rs (mul:rs cos-y zz) (mul:rs sin-y (add:rs (mul:rs sin-z yy) (mul:rs cos-z xx))))) (mul:rs cos-x (add:rs (mul:rs cos-z yy) (mul:rs sin-z xx))))
  =/  dz  (sub:rs (mul:rs cos-x (add:rs (mul:rs cos-y zz) (mul:rs sin-y (add:rs (mul:rs sin-z yy) (mul:rs cos-z xx))))) (mul:rs sin-x (add:rs (mul:rs cos-z yy) (mul:rs sin-z xx))))
  =/  bx  (add:rs (mul:rs ez (div:rs dx dz)) ex)
  =/  by  (add:rs (mul:rs ez (div:rs dy dz)) ey)
  ::~&  "{<ax>} {<ay>} {<az>} -> {<bx>} {<by>}"
  ~[bx by]
--
