::  Rössler Strange Attractor System
::  
::  By ~lagrev-nocfep (Neal Davis)
::
::  Solve the Rössler system using a Runge-Kutta fourth-order solution method.
::
::  dx/dt = -y - z
::  dy/dt = x + a*y
::  dz/dt = b + z*(x-c)
::
::  a = b = 0.1
::  c variable
::
|%
++  genxyzs
  ::  Produce a tang of the data from `genxyz`.
  |=  [nt=@ud dt=@rs]  ^-  tang
  %-  flop
  =/  x  .1
  =/  y  .1
  =/  z  .1
  =/  t  .0
  =/  index  1
  =/  acc=tang  ~
  |-  ^-  tang
  ?:  =(index +(nt))  acc
  =/  vals  (genxyz x y z index dt)
  $(index +(index), t -.vals, x +<.vals, y +>-.vals, z +>+<.vals, acc `tang`[leaf+"{<t>} {<x>} {<y>} {<z>}, " acc])
++  genxyz-pts
  ::  Produce a list of the data from `genxyz`.
  |=  [nt=@ud dt=@rs]  ^-  (list (list @rs))
  %-  flop
  =/  x  .1
  =/  y  .1
  =/  z  .1
  =/  t  .0
  =/  index  1
  =/  acc=(list (list @rs))  ~
  |-  ^-  (list (list @rs))
  ?:  =(index +(nt))  acc
  =/  vals  (genxyz x y z index dt)
  =/  tt  -.vals
  =/  xx  +<.vals
  =/  yy  +>-.vals
  =/  zz  +>+<.vals
  $(index +(index), t tt, x xx, y yy, z zz, acc `(list (list @rs))`[vals acc])
++  genxyz
  ::  Produce each time step as a list of `@rs`.
  |=  [x=@rs y=@rs z=@rs nt=@ud dt=@rs]  ^-  (list @rs)
  =/  tt  (mul:rs (sun:rs nt) dt)
  =/  a  .1e-1
  =/  b  .1e-1
  =/  c  .14
  =/  xx  (dxdt x y z tt a b c dt)
  =/  yy  (dydt x y z tt a b c dt)
  =/  zz  (dzdt x y z tt a b c dt)
  ~[tt xx yy zz]
::  For uniformity and future flexibility, all function signatures are the same.
++  xf
  ::  dx/dt = -y - z
  |=  [x=@rs y=@rs z=@rs t=@rs a=@rs b=@rs c=@rs]  ^-  @rs
  (sub:rs (sub:rs .0 y) z)
++  yf
  ::  dy/dt = x + a*y
  |=  [x=@rs y=@rs z=@rs t=@rs a=@rs b=@rs c=@rs]  ^-  @rs
  (add:rs x (mul:rs a y))
++  zf
  ::  dz/dt = b + z*(x-c)
  |=  [x=@rs y=@rs z=@rs t=@rs a=@rs b=@rs c=@rs]  ^-  @rs
  (add:rs b (mul:rs z (sub:rs x c)))
++  dxdt
  ::  RK4 implementation of dx/dt
  |=  [x=@rs y=@rs z=@rs t=@rs a=@rs b=@rs c=@rs dt=@rs]  ^-  @rs
  =/  k1  (xf x y z t a b c)
  =/  k2x  (add:rs x (mul:rs dt (div:rs k1 .2)))
  =/  k2y  (add:rs y (mul:rs dt (div:rs k1 .2)))
  =/  k2z  (add:rs z (mul:rs dt (div:rs k1 .2)))
  =/  k2t  (add:rs t (div:rs dt .2))
  =/  k2  (xf k2x k2y k2z k2t a b c)
  =/  k3x  (add:rs x (mul:rs dt (div:rs k2 .2)))
  =/  k3y  (add:rs y (mul:rs dt (div:rs k2 .2)))
  =/  k3z  (add:rs z (mul:rs dt (div:rs k2 .2)))
  =/  k3t  (add:rs t (div:rs dt .2))
  =/  k3  (xf k3x k3y k3z k3t a b c)
  =/  k4x  (add:rs x (mul:rs dt k3))
  =/  k4y  (add:rs y (mul:rs dt k3))
  =/  k4z  (add:rs z (mul:rs dt k3))
  =/  k4t  (add:rs t dt)
  =/  k4  (xf k4x k4y k4z k4t a b c)
  (add:rs x (div:rs (mul:rs dt (add:rs k1 (add:rs (mul:rs k2 .2) (add:rs (mul:rs k3 .2) k4)))) .6))
++  dydt
  ::  RK4 implementation of dy/dt
  |=  [x=@rs y=@rs z=@rs t=@rs a=@rs b=@rs c=@rs dt=@rs]  ^-  @rs
  =/  k1  (yf x y z t a b c)
  =/  k2x  (add:rs x (mul:rs dt (div:rs k1 .2)))
  =/  k2y  (add:rs y (mul:rs dt (div:rs k1 .2)))
  =/  k2z  (add:rs z (mul:rs dt (div:rs k1 .2)))
  =/  k2t  (add:rs t (div:rs dt .2))
  =/  k2  (yf k2x k2y k2z k2t a b c)
  =/  k3x  (add:rs x (mul:rs dt (div:rs k2 .2)))
  =/  k3y  (add:rs y (mul:rs dt (div:rs k2 .2)))
  =/  k3z  (add:rs z (mul:rs dt (div:rs k2 .2)))
  =/  k3t  (add:rs t (div:rs dt .2))
  =/  k3  (yf k3x k3y k3z k3t a b c)
  =/  k4x  (add:rs x (mul:rs dt k3))
  =/  k4y  (add:rs y (mul:rs dt k3))
  =/  k4z  (add:rs z (mul:rs dt k3))
  =/  k4t  (add:rs t dt)
  =/  k4  (yf k4x k4y k4z k4t a b c)
  (add:rs y (div:rs (mul:rs dt (add:rs k1 (add:rs (mul:rs k2 .2) (add:rs (mul:rs k3 .2) k4)))) .6))
++  dzdt
  ::  RK4 implementation of dz/dt
  |=  [x=@rs y=@rs z=@rs t=@rs a=@rs b=@rs c=@rs dt=@rs]  ^-  @rs
  =/  k1  (zf x y z t a b c)
  =/  k2x  (add:rs x (mul:rs dt (div:rs k1 .2)))
  =/  k2y  (add:rs y (mul:rs dt (div:rs k1 .2)))
  =/  k2z  (add:rs z (mul:rs dt (div:rs k1 .2)))
  =/  k2t  (add:rs t (div:rs dt .2))
  =/  k2  (zf k2x k2y k2z k2t a b c)
  =/  k3x  (add:rs x (mul:rs dt (div:rs k2 .2)))
  =/  k3y  (add:rs y (mul:rs dt (div:rs k2 .2)))
  =/  k3z  (add:rs z (mul:rs dt (div:rs k2 .2)))
  =/  k3t  (add:rs t (div:rs dt .2))
  =/  k3  (zf k3x k3y k3z k3t a b c)
  =/  k4x  (add:rs x (mul:rs dt k3))
  =/  k4y  (add:rs y (mul:rs dt k3))
  =/  k4z  (add:rs z (mul:rs dt k3))
  =/  k4t  (add:rs t dt)
  =/  k4  (zf k4x k4y k4z k4t a b c)
  (add:rs z (div:rs (mul:rs dt (add:rs k1 (add:rs (mul:rs k2 .2) (add:rs (mul:rs k3 .2) k4)))) .6))
--
