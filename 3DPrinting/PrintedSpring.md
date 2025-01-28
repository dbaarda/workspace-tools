

PLA properties;

roh = 1.25 Mg/m^3 AKA g/cm^2 == x water
E = 3.5GPa
smax = 35 MPa (some quote about 70 MPa)
G = 2.4 GPa
stress at yield 1~5% (2%?)
stress at break 10~40% (20%?)


I = w*h^3/12
S = w*h^2/6
M = P*L
d = P*L^3/(3*E*I) = P*L^3 / (3*E*w*h^3/12) = 4*P*L^3/(E*w*h^3)
s = P*L/S =P*L/(w*h^2/6) = 6*P*L/(w*h^2)

s = 3/2*E * h / L^2 * d

P = d*E*w*h^3/(4*L^3) = 1/4*E*w * (h/L)^3 * d


dmax = smax / (3/2*E * h / L^2) = smax * L^2 / (3/2*E * h)

Pmax = 1/4*E*w * (h/L)^3 * dmax
     = 1/4*E*w * (h/L)^3 * smax * L^2 / (3/2*E * h)
     = 1/4*w * h^3/L * smax / (3/2*h)
     = 1/6*w * smax * h^2/L
     
Pmax/smax = 1/6*w*h^2/L
  
dmax/Pmax = 4*L^3/(E*w*h^3) = 4/(E*w) * (L/h)^3

for a given displacement d0 and stress s0, how can we maximize P?


For s0:

s0 = P*6*L/(w*h^2)
P = s0*w*h^2/(6*L) = K1 * h^2/L

So for a fixed stress, the load increases with h^2/L.

For d0:

d0 = 4*P*L^3/(E*w*h^3)
P = d0*E*w*h^3/(4*L^3) = K2 * (h/L)^3

For a fixed displacement, the load increases with (h/L)^3


P as a function of d:

P = E*w*h^3/(4*L^3) * d = K3 * (h/L)^3 * d

s as a function of d

s = 3/2*E * h / L^2 * d = K4 * h/L^2 * d

So if you keep the ratio of h/L constant, the displacement for a given load
stays the same. However, if you keep h/L constant, the stress s decreases
linearly with L.

So h/L defines how much force/displacement you get, and then increasing h and
L together linearly decreases your stress. Both force and stress increase
linearly with displacement.


# Angle of bend

thinking of stress and displacement in terms of bend-angle of a straight beam. We have:

  L = R*a
  dL = h/2 * a
  R = L/a
  s = dL/L = h/2 * a / L
  d = (R+h/2)*(1-cos(a)) = (L/a + h/2)*(1-cos(a))
  
  a = 2*s*L/h
  
Where:

  L is length
  h is thickness
  a is bend-angle in radians
  R is radius of bend for the mid-chord
  dL is change in length at the inner/outter edge.
  s is max stress
  d is perpendicular displacement of the outer edge end.
  
Note for small a we can approximate with `cos(a) ~= 1 - a^2/2` to calculate `a`
from a given `L`, `h`, and `d`:
  
  1 - cos(a) ~= a^2/2
  d ~= (L/a + h/2)*a^2/2
     = L*a/2 + h*a^2/4
  4*d = h*a^2 + 2*L*a
  0 = h*a^2 + 2*L*a - 4*d
  
using solution to `0 = a*x^2 + b*x + c` of (-b +- sqrt(b^2 - 4*a*c))/(2*a)
  
  a = (-(2*L) +- sqrt((2*L)^2 + 4*h*(4*d)))/(2*h)
    = (-2*L +- sqrt(4*L^2 + 16*h*d))/(2*h)
    = (-L +- sqrt(L^2 + 4*h*d))/h

Example:

  L = 4.5mm
  h = 1.5mm
  a = 11.5deg = 0.2rad
  s = h/2 * a / L = 0.033 (3.3%)
  d = (L/a + h/2)*(1-cos(a)) ~= a/2*(L + h/2*a) = 0.47mm


