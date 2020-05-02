::  Transcendental functions library for Hoon, compatible with @rs
::
=/  tau  .6.28318530717
=/  pi   .3.14159265358
=/  e    .2.718281828
=/  rtol  .1e-5
|%
++  factorial
  |=  x=@rs  ^-  @rs
  =/  t=@rs  .1
  |-  ^-  @rs
  ?:  =(x .1)
    t
  $(x (sub:rs x .1), t (mul:rs t x))
++  absolute
  |=  x=@rs  ^-  @rs
  ?:  (gth:rs x .0)
    x
  (sub:rs .0 x)
++  exp
  |=  x=@rs  ^-  @rs
  =/  rtol  .1e-5
  =/  p   .1
  =/  po  .-1
  =/  i   .1
  |-  ^-  @rs
  ?:  (lth:rs (absolute (sub:rs po p)) rtol)
    p
  $(i (add:rs i .1), p (add:rs p (div:rs (pow-n x i) (factorial i))), po p)
++  pow-n
  ::  restricted pow, based on integers only
  |=  [x=@rs n=@rs]  ^-  @rs
  ?:  =(n .0)  .1
  =/  p  x
  |-  ^-  @rs
  ?:  (lth:rs n .2)
    p
  ::~&  [n p]
  $(n (sub:rs n .1), p (mul:rs p x))
++  sine
  ::  sin x = x - x^3/3! + x^5/5! - x^7/7! + x^9/9! - ...
  |=  x=@rs  ^-  @rs
  =/  rtol  .1e-5
  =/  p   .0
  =/  po  .-1
  =/  i   .0
  |-  ^-  @rs
  ?:  (lth:rs (absolute (sub:rs po p)) rtol)
    p
  =/  ii  (add:rs (mul:rs .2 i) .1)
  =/  term  (mul:rs (pow-n .-1 i) (div:rs (pow-n x ii) (factorial ii)))
  ::~&  [i ii term p]
  $(i (add:rs i .1), p (add:rs p term), po p)
++  cosine
  ::  cos x = 1 - x^2/2! + x^4/4! - x^6/6! + x^8/8! - ...
  |=  x=@rs  ^-  @rs
  =/  rtol  .1e-5
  =/  p   .1
  =/  po  .-1
  =/  i   .1
  |-  ^-  @rs
  ?:  (lth:rs (absolute (sub:rs po p)) rtol)
    p
  =/  ii  (mul:rs .2 i)
  =/  term  (mul:rs (pow-n .-1 i) (div:rs (pow-n x ii) (factorial ii)))
  ::~&  [i ii term p]
  $(i (add:rs i .1), p (add:rs p term), po p)
++  tangent
  ::  tan x = sin x / cos x
  |=  x=@rs  ^-  @rs
  (div:rs (sine x) (cosine x))
  ::  don't worry about domain errors right now, this is lazy trig
++  log-e-2
  ::  natural logarithm, only converges for z < 2
  |=  z=@rs  ^-  @rs
  =/  rtol  .1e-5
  =/  p   .0
  =/  po  .-1
  =/  i   .1
  |-  ^-  @rs
  ?:  (lth:rs (absolute (sub:rs po p)) rtol)
    p
  =/  ii  (add:rs .1 i)
  =/  term  (mul:rs (pow-n .-1 (add:rs .1 i)) (div:rs (pow-n (sub:rs z .1) i) i))
  ::~&  [i ii term p]
  $(i (add:rs i .1), p (add:rs p term), po p)
++  log-e
  ::  natural logarithm, z > 0
  |=  z=@rs  ^-  @rs
  =/  rtol  .1e-5
  =/  p   .0
  =/  po  .-1
  =/  i   .0
  |-  ^-  @rs
  ?:  (lth:rs (absolute (sub:rs po p)) rtol)
    (mul:rs (div:rs (mul:rs .2 (sub:rs z .1)) (add:rs z .1)) p)
  =/  term1  (div:rs .1 (add:rs .1 (mul:rs .2 i)))
  =/  term2  (mul:rs (sub:rs z .1) (sub:rs z .1))
  =/  term3  (mul:rs (add:rs z .1) (add:rs z .1))
  =/  term  (mul:rs term1 (pow-n (div:rs term2 term3) i))
  ::~&  [i term p]
  $(i (add:rs i .1), p (add:rs p term), po p)
++  pow
  ::  general power, based on logarithms
  ::  x^n = exp(n ln x)
  |=  [x=@rs n=@rs]  ^-  @rs
  (exp (mul:rs n (log-e x)))
--
