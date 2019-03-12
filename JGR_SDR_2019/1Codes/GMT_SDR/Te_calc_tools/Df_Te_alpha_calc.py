import math
def Df_Te():
    Df = float(input("enter the distance between axis and flat dike boundary/(steady SDRs geometry) in [km]: "))
    Df = Df * 1000; # convert the km into meter to fit into the units for the equation
    E = 5.0*10**10 #Pa
    g = 10 # m/s**2
    mu = 0.25
    delrhoC = float(input("enter delta rho C: (zero for pure lava, 1 for pure sediment, 2 for air) in [kg/m^3]:  "))

    if delrhoC==0:
        delrhoC = 200
    elif delrhoC==1:
        delrhoC = 600
    elif delrhoC ==2:
        delrhoC = 3000
        print("delrhoC: ", delrhoC)

    kappa = (E / (3 * delrhoC * g * (1 - mu**2)))**0.25 #[m]
    
    alpha = 2.0 * Df / math.pi  #[m]
    Te = (alpha / kappa)**(4.0/3.0) #[m]
    print ("Xf is:", Df, "[m]", "alpha_Xf is:", alpha, "[m]", "and Te_Xf is:", Te, "[m]", "kappa is:", kappa, "[m]") 
    return

Df_Te()
