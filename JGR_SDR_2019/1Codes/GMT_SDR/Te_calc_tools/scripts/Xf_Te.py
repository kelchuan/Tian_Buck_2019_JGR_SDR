# The Xf_Te() Afunction do two things:
# 1. it initializes parameters like: E, g, mu, delta_rho_c
# 2. it calculate alpha_xf and Te_xf from observation of Xf
import math
def Xf_Te():
    print("-------------------------------")
    print("begining of function Xf_Te()")
    print("-------------------------------")
    # 1. initialization/ parameter input
    E = 7.5*10**10 # Young's modulus [Pa]
    g = 10 # gravitational acc [m/s**2]
    mu = 0.25 # poisson's ratio
    delrhoC = 0; # type of delta_rho_C
    #delrhoC = float(input("enter delta rho C: (zero for pure lava, 1 for pure sediment, 2 for air) in [kg/m^3]:  "))
    if delrhoC==0:
        delrhoC = 200
    elif delrhoC==1:
        delrhoC = 600
    elif delrhoC ==2:
        delrhoC = 3000
    print("delrhoC: ", delrhoC)

    delrhoD = 200; # density difference between solidified and fluid dike
    
    # 2. Calculate alpha_xf and Te_xf from observable Xf
    
    Xf = float(input("enter the distance between axis and flat dike boundary/(steady SDRs geometry) in [km]: "))
    Xf = Xf * 1000; # convert the km into meter to fit into the units for the equation
    alpha_xf = 2.0 * Xf / math.pi  #[m]
    kappa = (E / (3 * delrhoC * g * (1 - mu**2)))**0.25 #coeficient between alpha_xf and Te_xf**0.75 [m] 
    Te_xf = (alpha_xf / kappa)**(4.0/3.0) #[m]

    coeff  = 2 * delrhoD / delrhoC**0.75 * (E / (3 * g * (1 - mu**2)))**-0.25
    print ("Xf is:", Xf, "[m]", "alpha_xf is:", alpha_xf, "[m]", "and Te_xf is:", Te_xf, "[m]", "kappa is:", kappa, "[m]")
    print("-------------------------------")
    print("end of function Xf_Te()")
    print("-------------------------------")
    return kappa,alpha_xf,Te_xf,coeff
