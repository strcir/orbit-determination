from __future__ import division 
from visual import *

NStars = 6
DecDeg = [37,37,37,37,37,37]
DecMin = [24,20,15,21,17,16]
DecSec = [52.85565,16.425,45.811,51.077,1.449,11.064]
RAHour = [16,15,15,15,16,15]
RAMin = [00,59,59,59,00,58]
RASec = [3.14889,41.1028,15.7558,17.2066,3.7434,21.3338]
Xpos =  [328.1441332,440.167790357,569.9277,564.2337,321.9286,851.7903]
Ypos = [159.8388145,281.338633762,401.0066,242.1154,364.1041,393.8594]
Xobjpos = 653
Yobjpos = 271
Dec = [0,0,0,0,0,0]
RA = [0,0,0,0,0,0]
Dec1 = [0,0,0,0,0,0]
RA1 = [0,0,0,0,0,0]
Xpos1 = [0,0,0,0,0,0]
Ypos1 = [0,0,0,0,0,0]
oo = [0,0,0,0,0,0]
pp = [0,0,0,0,0,0]
r1 = [0,0,0,0,0,0]
r2 = [0,0,0,0,0,0]

AvgDec = 0
AvgRA = 0
AvgX = 0
AvgY = 0
o = 0
p = 0
Starz1 = 0
Starz2 = 0
Objectz = 0
f = [0,0,0,0,0,0,0,0]
c1 = pi/180

#Convert RA, Dec to decimal, radians

for i in range(NStars):
    Dec[i] = DecDeg[i] + DecMin[i]/60 + DecSec[i]/3600
    RA[i] = RAHour [i] + RAMin[i]/60 + RASec[i]/3600
    Dec[i] = pi/180 * Dec[i]
    RA[i] = pi/12*RA[i]

for i in range(NStars): #calculate averages of Dec, RA, X and Y
    AvgDec = AvgDec + Dec[i]
    AvgRA = AvgRA + RA[i]
    AvgX = AvgX + Xpos[i]
    AvgY = AvgY + Ypos[i]

AvgDec = AvgDec/NStars
AvgRA = AvgRA/NStars
AvgX = (AvgX/NStars)/1200
AvgY = (AvgY/NStars)/1200

for i in range(NStars): #calculate differences for Least Squares
    RA1[i] = RA[i]
    Dec1[i] = Dec[i]
    Xpos1[i] = AvgX - Xpos[i]/1200
    Ypos1[i] = Ypos[i]/1200 - AvgY

X2 = AvgX - Xobjpos/1200.0
Y2 = Yobjpos/1200.0 - AvgY
D1 = sin(AvgDec)
D2 = cos(AvgDec)

for i in range(NStars):
    D3 = cos(Dec1[i])/sin(Dec1[i]) #convert RA and Dec to coordinates on idealized plate
    D4 = D1 + D2*D3*cos(RA1[i]-AvgRA)
    oo[i] = (D2 - D1*D3*cos(RA1[i]-AvgRA))/D4 - Ypos1[i]
    pp[i] = D3*sin(RA1[i]-AvgRA)/D4-Xpos1[i]
    o = o + oo[i]
    p = p + pp[i]
    f[0] = f[0] + Xpos1[i]*Xpos1[i]
    f[1] = f[1] + Ypos1[i]*Ypos1[i]
    f[3] = f[3] + Xpos1[i]*Ypos1[i]
    f[4] = f[4] + Xpos1[i]*pp[i]
    f[5] = f[5] + Ypos1[i]*pp[i]
    f[6] = f[6] + Xpos1[i]*oo[i]
    f[7] = f[7] + Ypos1[i]*oo[i]

D5 = (f[0]*f[1] - f[3]*f[3])
a0 = (f[4]*f[1] - f[3]*f[5])/D5
b0 = (f[0]*f[5] - f[4]*f[3])/D5
c0 = p / NStars
d0 = (f[6]*f[1]-f[3]*f[7])/D5
e0=(f[0]*f[7]-f[6]*f[3])/D5

f0 = o / NStars
t1 = D1 / D2
x1 = (a0 + 1)*X2 +b0*Y2 + c0
e1 = d0*X2 + (e0 + 1)*Y2 + f0
tt= atan((x1 / D2)/(1 -e1*t1)) + AvgRA
#Object Position
R1obj = tt/15.0/c1
R2obj = atan((e1 +t1)*cos(tt-AvgRA)/(1-e1*t1))/c1
wt=cos(c1*R2obj)*cos(c1*R2obj)

RAobjHour=int(R1obj)
RAobjMin= (R1obj-RAobjHour)*60
RAobjSec = (RAobjMin-int(RAobjMin))*60
RAobjMin=int(RAobjMin)
print "Object RA",RAobjHour,RAobjMin, RAobjSec
DecobjHour = int(R2obj)
DecobjMin = (R2obj-DecobjHour)*60
DecobjSec = (DecobjMin-int(DecobjMin))*60
DecobjMin = int(DecobjMin)
print "Object DEC", DecobjHour, abs(DecobjMin), abs(DecobjSec)

print "Star Residuals RA        DEC"
for i in range(NStars): #Generate Residuals for Stars
    x1 = (a0+1) * Xpos1[i] + b0*Ypos1[i]+c0
    e1 = d0*Xpos1[i]+(e0+1)*Ypos1[i]+f0
    tt=atan((x1/D2)/(1-e1*t1))+AvgRA
    r1[i]=(tt-RA1[i])/c1*3600.0
    r2[i]=(atan((e1+t1)*cos(tt-AvgRA)/(1-e1*t1))-Dec1[i])/c1*3600.0
    Starz1 +=r1[i]*r1[i]
    Starz2 +=r2[i]*r2[i]
    print" ",r1[i]," ",r2[i]

Starz1 = sqrt(Starz1 / (NStars-1))
Starz2 = sqrt(Starz2 / (NStars-1))
print "Mean Asteroid Error (RA/DEC)", Starz1, Starz2

Objectz=sqrt((Starz1)*(Starz1)*wt+(Starz2)*(Starz2))
print "Object Standard Deviation", Objectz
    


