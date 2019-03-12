import math
Dkm = float(input("enter the distance between max curvature and min curvature: "))
Dkm = Dkm * 1000 # convert to meter
E = 5.0*10**10 #Pa
g = 10 # m/s**2
mu = 0.25
delrhoC = float(input("enter delta rho C: (zero for pure lava, 1 for pure sediment, 2 for air)"))

if delrhoC==0:
    delrhoC = 200
elif delrhoC==1:
    delrhoC = 600
elif delrhoC ==2:
    delrhoC = 3000
print("delrhoC: ", delrhoC)

kappa = (E / (3 * delrhoC * g * (1 - mu**2)))**0.25

alpha = 4. / 3. * Dkm / math.pi
Te = (alpha / kappa)**(4.0/3.0)
print ("Dkm is:", Dkm, "alpha is : ", alpha, "and Te is:  ", Te, "kappa is: ", kappa) 

