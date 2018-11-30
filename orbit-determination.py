from __future__ import division
from math import *
from visual import *
from numpy import *

################################################################################
################################################################################
############## Orbit Determination and Ephemeris Generation Program ############
##############                   Strahinja Ciric                    ############
################################################################################
################################################################################


################################################################################
#Define required functions                                                     #
################################################################################

#define function to resolve angle ambiguities
def angleambig(sina, cosa):
    if sina >= 0 and cosa >= 0:
        a = asin(sina)
    elif sina > 0 and cosa < 0:
        a = acos(cosa)
    elif sina < 0 and cosa < 0:
        a = -2 * asin(sina) + pi
    elif sina < 0 and cosa > 0:
        a = asin(sina)
    if a < 0:
        a = a + 2 * pi
    return a

#Define function to convert hours to decimal hours
def hmshrs(hr, minutes, sec):
    if hr < 0:
        minutes = minutes*(-1)
        sec = sec*(-1)
    minutes = minutes + (sec / 60)
    hr = hr + (minutes / 60)
    if hr >= 0:
        hr = fmod(hr, 360)
    else:
        hr = fmod(hr, 360) + 360
    return hr

#Define function to convert degrees to decimal degrees
def dmsdeg(deg, minutes, sec):
    if deg < 0:
        minutes = minutes*(-1)
        sec = sec*(-1)
    minutes = minutes + (sec / 60)
    deg = deg + (minutes / 60)
    return deg

#define function to convert decimal degrees into radians
def degrad(deg):
    rad = deg / 180 * pi
    twopi = 2 * pi
    if rad >= 0:
        rad = fmod(rad, twopi)
    else:
        rad = fmod(rad, twopi) + twopi
    return rad

#Define function to convert hours into radians
def hrsrad(hr):
    rad = hr / 12 * pi
    twopi = 2 * pi
    if rad >= 0:
        rad = fmod(rad, twopi)
    else:
        rad = fmod(rad, twopi) + twopi
    return rad

#define function to convert radians into decimal degrees
def raddeg(rad):
    deg = rad * 180 / pi
    if deg >= 0:
        deg = fmod(deg, 360)
    else:
        deg = fmod(deg, 360) + 360
    return deg

#define function to convert radians into hours
def radhrs(rad):
    hr = rad * 12 / pi
    if hr >= 0:
        hr = fmod(hr, 24)
    else:
        hr = fmod(hr, 24) + 24
    return hr

################################################################################
#Determine orbital elements                                                    #
################################################################################

print """------------------------------------
        Orbit determination
------------------------------------

This program determines the orbital elements of asteroid
4055 Magellan.

In the following program outputs, the subscripts -1, 0, and 1
refer to the first, second, and third observations respectively.
rho is the Earth-Magellan vector, R is the Earth-Sol vector, and
r is the Sol-Magellan vector. a_1 and a_3 are coefficients
relating the three r vectors; they satisfy the equation

r_0 = a_1 * r_-1 + a_3 * r_1.

rho hat refers to the unit vector in the rho direction.
"""

#Create arrays for measurements
ra = zeros(3)
dec = zeros(3)
time = zeros(3)
rhohatx = zeros(3)
rhohatyi = zeros(3)
rhohaty = zeros(3)
rhohatz = zeros(3)
R = zeros(3)

#Input initial values. Note that the values of the variables in this
#section will be overwritten by those input from OrbitDet.txt if the
#file input code block is not put into comment form.
ra[0]= hrsrad(hmshrs(15.0, 59.0, 25.7599738348))
ra[1]= hrsrad(hmshrs(15.0, 58.0, 59.86))
ra[2]= hrsrad(hmshrs(16.0, 3.0, 10.05))
dec[0]= degrad(dmsdeg(37.0, 18.0, 20.5911033298))
dec[1]= degrad(dmsdeg(35.0, 5.0, 31.8))
dec[2]= degrad(dmsdeg(31.0, 33.0, 55.4))
time[0]= 2455377.69499
time[1]= 2455386.64517
time[2]= 2455396.68779
for i in range(0, 3):
    rhohatx[i] = cos(ra[i])*cos(dec[i])
    rhohatyi[i] = sin(ra[i])*cos(dec[i])
    rhohatz[i] = sin(dec[i])

R1n = vector(-1.443910920559232E-01,  1.006302121603611E+00,
             -2.137449447053117E-05)
R0 = vector(-2.921341251003500E-01,  9.738073410270751E-01,
            -2.053008927441118E-05)
R1 = vector(-4.499398431342010E-01,  9.111913267305856E-01,
            -2.253067405596743E-05)

#Define constants
epsilon = degrad(23.439281)

#The following code block, when uncommented, will enable file input from
#OrbitDet.txt.
#Note that the results are for more accurate when this file input system
#is used.

data = open("OrbitDet.txt", "r")
time[0] = float(data.readline())
obs1n = data.readline()
obs1n = obs1n.split()
ra[0] = hrsrad(hmshrs(float(obs1n[0]), float(obs1n[1]), float(obs1n[2])))
dec[0] = degrad(dmsdeg(float(obs1n[3]), float(obs1n[4]), float(obs1n[5])))
R1n = data.readline()
R1n = R1n.split()
R1n = [float(R1n[0]), float(R1n[1]), float(R1n[2])]
R1ny = R1n[1]*cos(epsilon) + R1n[2]*sin(epsilon)
R1n[2] = R1n[1]*(-sin(epsilon)) + (R1n[2]*cos(epsilon))
R1n[1] = R1ny
R1n = vector(R1n[0], R1n[1], R1n[2])
linebreak = data.readline()
time[1] = float(data.readline())
obs0 = data.readline()
obs0 = obs0.split()
ra[1] = hrsrad(hmshrs(float(obs0[0]), float(obs0[1]), float(obs0[2])))
dec[1] = degrad(dmsdeg(float(obs0[3]), float(obs0[4]), float(obs0[5])))
R0 = data.readline()
R0 = R0.split()
R0 = [float(R0[0]), float(R0[1]), float(R0[2])]
R0y = R0[1]*cos(epsilon) + R0[2]*sin(epsilon)
R0[2] = R0[1]*(-sin(epsilon)) + (R0[2]*cos(epsilon))
R0[1] = R0y
R0 = vector(R0[0], R0[1], R0[2])
linebreak = data.readline()
time[2] = float(data.readline())
obs1 = data.readline()
obs1 = obs1.split()
ra[2] = hrsrad(hmshrs(float(obs1[0]), float(obs1[1]), float(obs1[2])))
dec[2] = degrad(dmsdeg(float(obs1[3]), float(obs1[4]), float(obs1[5])))
R1 = data.readline()
R1 = R1.split()
R1 = [float(R1[0]), float(R1[1]), float(R1[2])]
R1y = R1[1]*cos(epsilon) + R1[2]*sin(epsilon)
R1[2] = R1[1]*(-sin(epsilon)) + (R1[2]*cos(epsilon))
R1[1] = R1y
R1 = vector(R1[0], R1[1], R1[2])

#Change rho hat vector components to ecliptic coordinates
rhohaty[1] = rhohatyi[1]*cos(epsilon) + rhohatz[1]*sin(epsilon)
rhohaty[2] = rhohatyi[2]*cos(epsilon) + rhohatz[2]*sin(epsilon)
rhohaty[0] = rhohatyi[0]*cos(epsilon) + rhohatz[0]*sin(epsilon)
rhohatz[1] = rhohatyi[1]*(-sin(epsilon)) + (rhohatz[1]*cos(epsilon))
rhohatz[2] = rhohatyi[2]*(-sin(epsilon)) + (rhohatz[2]*cos(epsilon))
rhohatz[0] = rhohatyi[0]*(-sin(epsilon)) + (rhohatz[0]*cos(epsilon))

#Define rho hat vectors
rhohat0 = vector(rhohatx[1], rhohaty[1], rhohatz[1])
rhohat1 = vector(rhohatx[2], rhohaty[2], rhohatz[2])
rhohat1n = vector(rhohatx[0], rhohaty[0], rhohatz[0])

#Find k
k = 0.01720209894

#Find tau
tau1n = k*(time[0]-time[1])
tau0 = k*(time[2]-time[0])
tau1 = k*(time[2]-time[1])

#Find initial coefficients of r1n and r1
a1 = tau1/tau0
a3 = -tau1n/tau0

#Find initial rho
rho1n = ((a1*dot(cross(R1n, rhohat0), rhohat1)-
          dot(cross(R0, rhohat0), rhohat1)+
          a3*dot(cross(R1, rhohat0), rhohat1))/
         (a1*dot(cross(rhohat1n, rhohat0), rhohat1)))
rho0 = ((a1*dot(cross(rhohat1n, R1n), rhohat1)-
          dot(cross(rhohat1n, R0), rhohat1)+
          a3*dot(cross(rhohat1n, R1), rhohat1))/
         (-dot(cross(rhohat1n, rhohat0), rhohat1)))
rho1 = ((a1*dot(cross(rhohat0, R1n), rhohat1n)-
          dot(cross(rhohat0, R0), rhohat1n)+
          a3*dot(cross(rhohat0, R1), rhohat1n))/
         (a3*dot(cross(rhohat0, rhohat1), rhohat1n)))

#Find initial r
r1n = rho1n*rhohat1n-R1n
r1 = rho1*rhohat1-R1
r0 = rho0*rhohat0-R0

#Find the initial time derivative of r
rdot = ((r1-r0)/(tau1)-(r0-r1n)/(tau1n))/2

#Iterate f and g series
a1a = 0
a3a = 0
count = 0
print """Initial approximated values:
"""
print "a_1 =", a1
print "a_3 =", a3
print "r_0 =", r0
print "dr_0/dt =", rdot
print """
--iterated series to improve approximation--

Using f- and g-series, these values are more and more closely
determined until they are found to converge so that the
difference between a1 in successive loops is smaller than
0.0000000001. Values at each iteration are displayed below
in the following format:

( a1, previous loop's a1 ) ( a3, previous loop's a3 )
r0 vector
"""
while (abs(a1-a1a) > 0.0000000001):
    #Set standards for comparison for test to exit loop
    a1a = a1
    a3a = a3
    
    #Determine f and g
    f1n = (1-(tau1n**2/(2*mag(r0)**3))+
           (((tau1n**3)*dot(r0, rdot))/(2*mag(r0)**5))+
           ((tau1n**4)/24*(3/(mag(r0)**3)*(dot(rdot, rdot)/(mag(r0)**2)-1/
                                           mag(r0)**3)-15/(mag(r0)**7)
                           *(dot(r0, rdot))**2+1/(mag(r0)**6))))

    f1 = (1-tau1**2/(2*mag(r0)**3)+
           ((tau1**3)*(dot(r0, rdot)/(2*mag(r0)**5)))+
           ((tau1**4)/24*(3/(mag(r0)**3)*(dot(rdot, rdot)/(mag(r0)**2)-1/
                                           mag(r0)**3)-15/(mag(r0)**7)
                           *(dot(r0, rdot))**2+1/(mag(r0)**6))))
    g1n = (tau1n-((tau1n**3)/(6*mag(r0)**3))
           +((tau1n**4)*(dot(r0, rdot)))/(4*mag(r0)**5))
    g1 = (tau1-((tau1**3)/(6*mag(r0)**3))
          +((tau1**4)*(dot(r0, rdot)))/(4*mag(r0)**5))

    #Refine measurement of a
    a1 = g1/(f1n*g1-f1*g1n)
    a3 = -g1n/(f1n*g1-f1*g1n)

    #Refine measurements of rho nought and r
    rho1n = ((a1*dot(cross(R1n, rhohat0), rhohat1)-
              dot(cross(R0, rhohat0), rhohat1)+
              a3*dot(cross(R1, rhohat0), rhohat1))/
             (a1*dot(cross(rhohat1n, rhohat0), rhohat1)))
    rho0 = ((a1*dot(cross(rhohat1n, R1n), rhohat1)-
              dot(cross(rhohat1n, R0), rhohat1)+
              a3*dot(cross(rhohat1n, R1), rhohat1))/
             (-dot(cross(rhohat1n, rhohat0), rhohat1)))
    rho1 = ((a1*dot(cross(rhohat0, R1n), rhohat1n)-
              dot(cross(rhohat0, R0), rhohat1n)+
              a3*dot(cross(rhohat0, R1), rhohat1n))/
             (a3*dot(cross(rhohat0, rhohat1), rhohat1n)))
    r0 = rho0*rhohat0-R0
    r1n = rho1n*rhohat1n-R1n
    r1 = rho1*rhohat1-R1

    #Refine derviative measurement
    rdot = (r1-f1*r0)/(g1)
    print "(", a1, a1a, ")", "(", a3, a3a, ")"
    print "r0 =", r0

    #Increase iteration count
    count=count+1

#Print final values
print """
--end of iterated series--
"""
print"The final values determined after", count, "iterations are as"
print"""follows:
"""
print "rho_-1 =", rho1n*rhohat1n
print "rho_0 =", rho0*rhohat0
print "rho_1 =", rho1*rhohat1
print "r_-1 =", r1n
print "r_0 =", r0
print "r_1 =", r1
print "dr_0/dt =", rdot
print """
The orbital elements as determined from the given right
ascensions, declinations, and Earth-Sol vectors are as follows:
"""

#Find the semimajor axis
a = 1/(2/mag(r0)-(dot(rdot, rdot)))
print "Semimajor axis =", a, "AU"

#Find the eccentricity
e = sqrt(1-((mag(cross(r0, rdot)))**2)/(a))
print "Eccentricity =", e

#Find the inclination
h = cross(r0, rdot)
i = atan(sqrt((h[0])**2+(h[1])**2)/(h[2]))
print "Inclination =", raddeg(i), "degrees"

#Find the longitude of the ascending node
sinomegacap = h[0]/(mag(h)*sin(i))
cosomegacap = -h[1]/(mag(h)*sin(i))
omegacap = angleambig(sinomegacap, cosomegacap)
print "Longitude of ascending nodes =", raddeg(omegacap), "degrees"

#Find the argument of the perihelion
cosU = (r0[0]*cos(omegacap)+r0[1]*sin(omegacap))/(mag(r0))
sinU = r0[2]/((mag(r0))*sin(i))
U = angleambig(sinU, cosU)
cosnu = (((a*(1-e**2))/(mag(r0)))-1)/e
sinnu = ((a*(1-e**2)*(dot(r0, rdot)))/(mag(r0)*mag(h)))/e
nu = angleambig(sinnu, cosnu)
omegasmall = U - nu
print "Argument of perihelion =", raddeg(omegasmall), "degrees"

#Find the mean anomaly
E = acos(1/e*(1-mag(r0)/a))
M = E - e*sin(E)
if time[1]>=2455400.5:
    MtE = k/(a**(3/2))*(2455400.5-time[1])+M
if time[1]<2455400.5:
    MtE = k/(a**(3/2))*(2455400.5-time[1])-M
print "Mean anomaly at observation =", raddeg(M)
print "Reference mean anomaly =", raddeg(MtE)

################################################################################
#Generate an ephemeris                                                         #
################################################################################

print """
------------------------------------
        Ephemeris generation
------------------------------------
"""

#Set ephemeris time and solar vector
timeeph = 2455387.76666667
Reph = vector(-3.102607950592273E-01, 9.681870001238805E-01,
              -5.681569399685279E-05)

#Convert the solar vector to ecliptic coordinates
Rephy = Reph.y*cos(epsilon) - Reph.z*sin(epsilon)
Reph.z = Reph.y*(sin(epsilon)) + (Reph.z*cos(epsilon))
Reph.y = Rephy

#Find the mean anomaly at the ephemeris time
Meph = (k/(a**(3/2)))*(timeeph-2455400.5)+MtE

#Use Newton's method to find the eccentric anomaly
Eeph = Meph
Meph1 = 0.0
while (abs(Meph1-Meph) > 0.0000000001):
    Eeph = Eeph-((Meph-(Eeph-e*sin(Eeph)))/(e*cos(Eeph)-1.0))
    Meph1 = Eeph-e*sin(Eeph)

#Find the position vector in Cartesian coordinates
posCart = vector(a*cos(Eeph)-a*e, a*sqrt(1-e**2)*sin(Eeph), 0.0)

#Convert r to ecliptic coordinates
matrix111=cos(omegacap)
matrix112=cos(i)*-sin(omegacap)
matrix113=-sin(i)*-sin(omegacap)
matrix121=sin(omegacap)
matrix122=cos(i)*cos(omegacap)
matrix123=-sin(i)*cos(omegacap)
matrix131=0
matrix132=sin(i)
matrix133=cos(i)
matrix211=matrix111*cos(omegasmall)+matrix112*sin(omegasmall)
matrix212=matrix111*-sin(omegasmall)+matrix112*cos(omegasmall)
matrix213=matrix113
matrix221=matrix121*cos(omegasmall)+matrix122*sin(omegasmall)
matrix222=matrix121*-sin(omegasmall)+matrix122*cos(omegasmall)
matrix223=matrix123
matrix231=matrix131*cos(omegasmall)+matrix132*sin(omegasmall)
matrix232=matrix131*-sin(omegasmall)+matrix132*cos(omegasmall)
matrix233=matrix133
reclx=matrix211*posCart.x+matrix212*posCart.y+matrix213*posCart.z
recly=matrix221*posCart.x+matrix222*posCart.y+matrix223*posCart.z
reclz=matrix231*posCart.x+matrix232*posCart.y+matrix233*posCart.z
recleph=vector(reclx, recly, reclz)

#Convert r to equatorial coordinates
requx=reclx
requy=recly*cos(epsilon)+reclz*-sin(epsilon)
requz=recly*sin(epsilon)+reclz*cos(epsilon)
reph=vector(requx, requy, requz)

print "r =", reph

#Find the rho vector
rhoeph = reph + Reph

print "rho =", rhoeph

#Find the unit vector in thr rho direction
rhohateph = rhoeph/(mag(rhoeph))

print "rho hat =", rhohateph

#Find right ascension and declination
deceph = asin(rhohateph.z)
cosraeph = rhohateph.x/cos(deceph)
sinraeph = rhohateph.y/cos(deceph)
raeph = angleambig(sinraeph, cosraeph)
decephdeg = raddeg(deceph)
raephhrs = radhrs(raeph)
ramineph = (raephhrs-int(raephhrs))*60
raseceph = (ramineph-int(ramineph))*60
ramineph = int(ramineph)
raephhrs = int(raephhrs)
decmineph = (decephdeg - int(decephdeg))*60
decseceph = (decmineph - int(decmineph))*60
decmineph = int(decmineph)
decephdeg = int(decephdeg)

print "Declination =", decephdeg, decmineph, decseceph
print "Right Ascension =", raephhrs, ramineph, raseceph
