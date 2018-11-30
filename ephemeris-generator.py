from visual import *
from math import *
#from scipy import *
from Numeric import *

# functions
# min sec conversion
def convert(x, y, z):
    mintodeg = y / 60.0
    sectodeg = z / 3600.0
    minplussec = mintodeg + sectodeg
    if x < 0:
        x = x - minplussec
    else:
        x = x + minplussec
    return x

# f series
def f(tau, r0, v0):
    term1 = 1
    term2 = tau ** 2 / (2 * mag(r0) ** 3)
    term3 = tau ** 3 * dot(r0, v0) / (2 * mag(r0) ** 5)
    term4_1 = tau ** 4 / 24
    term4_2 = (3 / mag(r0) ** 3) * (dot(v0, v0) / mag(r0) ** 2 - 1 / mag(r0) ** 3)
    term4_3 = (15 / mag(r0) ** 7) * dot(v0, r0) ** 2
    term4_4 = 1 / mag(r0) ** 6
    value = term1 - term2 + term3 + term4_1 * (term4_2 - term4_3 + term4_4)
    return value

# g series
def g(tau, r0, v0):
    term1 = tau
    term2 = tau ** 3 / (6 * mag(r0) ** 3)
    term3 = tau ** 4 * dot(v0, r0) / (4 * mag(r0) ** 5)
    value = term1 - term2 + term3
    return value

# angle ambiguity
def angleambig(sin, cos):
    # Quadrant I
    if sin >= 0 and cos >= 0:
        angle = asin(sin)
    # Quadrant II
    elif sin >= 0 and cos <= 0:
        angle = pi - asin(sin)
    # Quadrant III
    elif sin <= 0 and cos <= 0:
        angle = pi - asin(sin)
    # Quadrant IV
    elif sin <= 0 and cos >= 0:
        angle =  2.0 * pi + asin(sin)
    return angle

# Get data
data = open("TrueODCalvin.txt", "r")
JD1 = data.readline()
JD1 = float(JD1)
coords1 = data.readline()
coords1 = coords1.split()
useless = data.readline()
JD2 = data.readline()
JD2 = float(JD2)
coords2 = data.readline()
coords2 = coords2.split()
useless = data.readline()
JD3 = data.readline()
JD3 = float(JD3)
coords3 = data.readline()
coords3 = coords3.split()

solarminusone = vector(-0.202887, 0.995827, 0)
solarzero = vector(-0.254176, 0.984430, 0)
solarplusone = vector(-0.351321, 0.953996, 0)

# RA's & Dec's
for x in range(0, 6):
    coords1[x] = float(coords1[x])
    coords2[x] = float(coords2[x])
    coords3[x] = float(coords3[x])
    
# Convert to radians & get rho hats (in equatorial)
RA1 = convert(coords1[0], coords1[1], coords1[2]) * pi / 12.0
RA2 = convert(coords2[0], coords2[1], coords2[2]) * pi / 12.0
RA3 = convert(coords3[0], coords3[1], coords3[2]) * pi / 12.0
Dec1 = convert(coords1[3], coords1[4], coords1[5]) * pi / 180.0
Dec2 = convert(coords2[3], coords2[4], coords2[5]) * pi / 180.0
Dec3 = convert(coords3[3], coords3[4], coords3[5]) * pi / 180.0
print RA1
rhohatminusone = vector(cos(RA1) * cos(Dec1), sin(RA1) * cos(Dec1), sin(Dec1))
rhohatzero = vector(cos(RA2) * cos(Dec2), sin(RA2) * cos(Dec2), sin(Dec2))
rhohatplusone = vector(cos(RA3) * cos(Dec3), sin(RA3) * cos(Dec3), sin(Dec3))

# Convert rho hats to ecliptic
ecliptictilt = 23.439281 * pi / 180.0

store1 = rhohatminusone.y
store2 = rhohatminusone.z
store3 = rhohatzero.y
store4 = rhohatzero.z
store5 = rhohatplusone.y
store6 = rhohatplusone.z

rhohatminusone.y = store1 * cos(ecliptictilt) + store2 * sin(ecliptictilt)
rhohatminusone.z = store2 * cos(ecliptictilt) - store1 * sin(ecliptictilt)

rhohatzero.y = store3 * cos(ecliptictilt) + store4 * sin(ecliptictilt)
rhohatzero.z = store4 * cos(ecliptictilt) - store3 * sin(ecliptictilt)

rhohatplusone.y = store5 * cos(ecliptictilt) + store6 * sin(ecliptictilt)
rhohatplusone.z = store6 * cos(ecliptictilt) - store5 * sin(ecliptictilt)

# Gaussian units
k = 0.01720209895
tauminusone = k * (JD1 - JD2)
tauzero = k * (JD3 - JD1)
tauplusone = k * (JD3 - JD2)

# r0 components in terms of taus
comp1 = tauplusone / tauzero
comp3 = -1 * tauminusone / tauzero

# rho vector magnitudes (stupidlongtyping)
minusonecross1 = cross(solarminusone, rhohatzero)
minusonecross2 = cross(solarzero, rhohatzero)
minusonecross3 = cross(solarplusone, rhohatzero)
minusonecross4 = cross(rhohatminusone, rhohatzero)
minusonedot1 = dot(minusonecross1, rhohatplusone)
minusonedot2 = dot(minusonecross2, rhohatplusone)
minusonedot3 = dot(minusonecross3, rhohatplusone)
minusonedot4 = dot(minusonecross4, rhohatplusone)
rhomagminusone = (comp1 * minusonedot1 - minusonedot2 + comp3 * minusonedot3) / (comp1 * minusonedot4)
rhomagminusone = abs(rhomagminusone)

zerocross1 = cross(rhohatminusone, solarminusone)
zerocross2 = cross(rhohatminusone, solarzero)
zerocross3 = cross(rhohatminusone, solarplusone)
zerocross4 = cross(rhohatminusone, rhohatzero)
zerodot1 = dot(zerocross1, rhohatplusone)
zerodot2 = dot(zerocross2, rhohatplusone)
zerodot3 = dot(zerocross3, rhohatplusone)
zerodot4 = dot(zerocross4, rhohatplusone)
rhomagzero = (comp1 * zerodot1 - zerodot2 + comp3 * zerodot3) / (-1 * zerodot4)
rhomagzero = abs(rhomagzero)

plusonecross1 = cross(rhohatzero, solarminusone)
plusonecross2 = cross(rhohatzero, solarzero)
plusonecross3 = cross(rhohatzero, solarplusone)
plusonecross4 = cross(rhohatzero, rhohatplusone)
plusonedot1 = dot(plusonecross1, rhohatminusone)
plusonedot2 = dot(plusonecross2, rhohatminusone)
plusonedot3 = dot(plusonecross3, rhohatminusone)
plusonedot4 = dot(plusonecross4, rhohatminusone)
rhomagplusone = (comp1 * plusonedot1 - plusonedot2 + comp3 * plusonedot3) / (comp3 * plusonedot4)
rhomagplusone = abs(rhomagplusone)
print "rhomagplusone is:", rhomagplusone
# rho vectors
rhominusone = rhomagminusone * rhohatminusone
rhozero = rhomagzero * rhohatzero
rhoplusone = rhomagplusone * rhohatplusone
print "rhominusone is:", rhominusone, "rhozero is:", rhozero, "rhoplusone is:", rhoplusone
# first r vectors
rminusone = rhominusone - solarminusone
rplusone = rhoplusone - solarplusone
rzero = rhozero - solarzero
print "rminusone is:", rminusone, "rplusone:", rplusone, "rzero is:", rzero

# first rdot vector
rdotzero = ((rplusone - rzero) / tauplusone + (rzero - rminusone) / (-1 * tauminusone)) / 2
print rdotzero, "RDotZero"
# loop 'o doom
check1 = 0
check2 = 0
count = 0
print comp1, comp3
while abs(check1 - comp1) > 1e-14 or abs(check3 - comp3) > 1e-14:
    check1 = comp1
    check3 = comp3
    compdenom = f(tauminusone, rzero, rdotzero) * g(tauplusone, rzero, rdotzero) - f(tauplusone, rzero, rdotzero) * g(tauminusone, rzero, rdotzero)
    comp1 = g(tauplusone, rzero, rdotzero) / compdenom 
    comp3 = -1 * g(tauminusone, rzero, rdotzero) / compdenom
    rdotzero = (f(tauplusone, rzero, rdotzero) * rminusone - f(tauminusone, rzero, rdotzero) * rplusone) / (f(tauplusone, rzero, rdotzero) * g(tauminusone, rzero, rdotzero) - f(tauminusone, rzero, rdotzero) * g(tauplusone, rzero, rdotzero))

    print comp1, comp3
    
    #print rhomagminusone, rhomagzero, rhomagplusone
    
    minusonecross1 = cross(solarminusone, rhohatzero)
    minusonecross2 = cross(solarzero, rhohatzero)
    minusonecross3 = cross(solarplusone, rhohatzero)
    minusonecross4 = cross(rhohatminusone, rhohatzero)
    minusonedot1 = dot(minusonecross1, rhohatplusone)
    minusonedot2 = dot(minusonecross2, rhohatplusone)
    minusonedot3 = dot(minusonecross3, rhohatplusone)
    minusonedot4 = dot(minusonecross4, rhohatplusone)
    rhomagminusone = (comp1 * minusonedot1 - minusonedot2 + comp3 * minusonedot3) / (comp1 * minusonedot4)
    rhomagminusone = abs(rhomagminusone)
    
    zerocross1 = cross(rhohatminusone, solarminusone)
    zerocross2 = cross(rhohatminusone, solarzero)
    zerocross3 = cross(rhohatminusone, solarplusone)
    zerocross4 = cross(rhohatminusone, rhohatzero)
    zerodot1 = dot(zerocross1, rhohatplusone)
    zerodot2 = dot(zerocross2, rhohatplusone)
    zerodot3 = dot(zerocross3, rhohatplusone)
    zerodot4 = dot(zerocross4, rhohatplusone)
    rhomagzero = (comp1 * zerodot1 - zerodot2 + comp3 * zerodot3) / (-1 * zerodot4)
    rhomagzero = abs(rhomagzero)
    
    plusonecross1 = cross(rhohatzero, solarminusone)
    plusonecross2 = cross(rhohatzero, solarzero)
    plusonecross3 = cross(rhohatzero, solarplusone)
    plusonecross4 = cross(rhohatzero, rhohatplusone)
    plusonedot1 = dot(plusonecross1, rhohatminusone)
    plusonedot2 = dot(plusonecross2, rhohatminusone)
    plusonedot3 = dot(plusonecross3, rhohatminusone)
    plusonedot4 = dot(plusonecross4, rhohatminusone)
    rhomagplusone = (comp1 * plusonedot1 - plusonedot2 + comp3 * plusonedot3) / (comp3 * plusonedot4)
    rhomagplusone = abs(rhomagplusone)

    rhominusone = rhomagminusone * rhohatminusone
    rhoplusone = rhomagplusone * rhohatplusone
    rhozero = rhomagzero * rhohatzero

    rminusone = rhominusone - solarminusone
    rplusone = rhoplusone - solarplusone

    rzero = comp1 *  rminusone + comp3 * rplusone

    count = count + 1
print "--------------------------------------------------------"
print count, "iterations"
print comp1, "a1"
print comp3, "a3"
print rzero, "r[0]"
print rdotzero, "rDot[0]"
print rhozero, "Rho[0]"

# The orbital elements

h = cross(rzero, rdotzero)

# Semimajor Axis

a = 1 / (2 / mag(rzero) - dot(rdotzero, rdotzero))

# Eccentricity

e = sqrt(1 - mag(h) ** 2 / (a))

# Inclination

i = atan(((h.x ** 2 + h.y ** 2) ** (0.5)) / h.z)
displayi = i * 180 / pi

# Longitude of the Ascending Node

sinom = h.x / (mag(h) * sin(i))
cosom = - h.y / (mag(h) * sin(i))
bigomega = angleambig(sinom, cosom)
displaybigomega = bigomega * 180 / pi

# Argument of Perihelion

cosu = (rzero.x * cos(bigomega) + rzero.y * sin(bigomega)) / mag(rzero)
sinu = rzero.z / (mag(rzero) * sin(i))
ecosnu = (a * (1 - e ** 2) / mag(rzero) - 1)
esinnu = a * (1 - e ** 2) * dot(rzero, rdotzero) / (mag(h) * mag(rzero))
cosnu = ecosnu / e
sinnu = esinnu / e
u = angleambig(sinu, cosu)
nu = angleambig(sinnu, cosnu)
littleomega = u - nu
if littleomega < 0:
    littleomega = 2 * pi + littleomega
displaylittleomega = littleomega * 180 / pi

# Mean Anomaly

n = k / (a ** 1.5)
diff = 2455400.5 - JD2
print JD2, "JD2"
E = acos((1 / e) * (1 - mag(rzero) / a))
M = E - e * sin(E)
bah = M * 180 / pi
if JD2 > 2455400.5:
    manom = n * diff + M
if JD2 < 2455400.5:
    manom = n * diff - M
while manom < 0:
    manom = manom + 2 * pi
while manom > 2 * pi:
    manom = manom - 2 * pi
manom = manom * 180 / pi
print "--------------------------------------------------------"
print "Elements:"
print "--------------------------------------------------------"
print "Semimajor Axis:", a
print "Eccentricity:", e
print "Inclination:", displayi
print "Longitude of the Ascending Node:", displaybigomega
print "Argument of Perihelion:", displaylittleomega
print "Mean Anomaly:", manom

#a  =  0.882859790401
#e  =  0.633549014396
#i  =  17.0851805457 * pi / 180
#bigomega  =  130.244179669 * pi / 180
#littleomega  =  309.20047405 * pi / 180
#manom  =  281.069174929

# Ephemeris generation

newM = 0.0
print "--------------------------------------------------------"
t = input("Input your fourth observation's JD: ")

# Earth-Sun Vector
# Constants
G = 6.67828e-11 / ((1.49597870692e11)**3)
K = ((3.154e7) ** 2)
# Earth Stats
earth_e = 0.016710220
earth_sma = 1.00000261
earth_period = sqrt(earth_sma ** 3 * K)
theta = 15.0180555556 * 86400.0 / earth_period * 2 * pi
earth_mass = 5.9742e24
earth_aphelion = earth_sma * (1 + earth_e)
earth_pos = vector(earth_aphelion * sin(theta), -earth_aphelion * cos(theta), 0)
earth_smina = earth_sma * sqrt(1 - earth_e ** 2)
earth_sweep = pi * earth_sma * earth_smina / earth_period
earth_initvel = 2 * earth_sweep / earth_aphelion
earth_vel = vector(earth_initvel * cos(theta), earth_initvel * sin(theta), 0)
# Sun Stats
sun_mass = 1.988435e30
# Init JD
JD = 2455384.0000000
# dt
dt = 33.0
# conditions for asigning solar vector values
solarfound = false
# find the vectors
while 1==1:
    earth_pos = earth_pos + earth_vel * dt
    earth_accel = -1 * (G * sun_mass / mag(earth_pos) ** 3) * earth_pos
    earth_vel = earth_vel + earth_accel * dt
    JD = JD + dt / 86400.0
    if JD >= t and solarfound == false:
        solar = -earth_pos
        solarfound = true
        break
n = k / a ** 1.5
solar = vector(solar.x, solar.y * cos(ecliptictilt) - solar.z * sin(ecliptictilt), solar.y * sin(ecliptictilt) + solar.z * cos(ecliptictilt))
origM = n * (t - 2455400.5) + manom * pi / 180

newtonx = origM
print newtonx
while abs(newM - origM) > 1e-9:
    newtonx = newtonx - ((origM - (newtonx - e * sin(newtonx))) / (e * cos(newtonx) - 1.0))
    newM = newtonx - e * sin(newtonx)
    print E, newtonx, newM, origM
print a, e, newtonx, "a,e,Newton"
posvec = mat([[a * cos(newtonx) - a * e], [a * sqrt(1 - e ** 2) * sin(newtonx)], [0]])
A = mat([[cos(bigomega), -sin(bigomega), 0], [sin(bigomega), cos(bigomega), 0], [0, 0, 1]])
B = mat([[1, 0, 0], [0, cos(i), -sin(i)], [0, sin(i), cos(i)]])
C = mat([[cos(littleomega), -sin(littleomega), 0], [sin(littleomega), cos(littleomega), 0], [0, 0, 1]])

posvec = C * posvec
posvec = B * posvec
posvec = A * posvec

D = mat([[1, 0, 0], [0, cos(ecliptictilt), -sin(ecliptictilt)], [0, sin(ecliptictilt), cos(ecliptictilt)]])

posvec = D * posvec
posvec = vector(take(posvec, [0])[0][0], take(posvec, [1])[0][0], take(posvec, [2])[0][0])
ephemrho = posvec + solar
ephemrhohat = ephemrho / mag(ephemrho)
dec = asin(ephemrhohat.z)
cosra = ephemrhohat.x / cos(dec)
sinra = ephemrhohat.y / cos(dec)
ra = angleambig(sinra, cosra)
ra = ra * 12 / pi
dec = dec * 180 / pi

print ra, dec
