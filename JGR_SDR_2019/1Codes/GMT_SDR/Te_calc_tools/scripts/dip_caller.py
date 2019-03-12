import math
import true_dip_with_vertical_exaggeration as TD

# Te estimation from ratio calculation with real depth
#def RD():
#    Ws = float(input("enter the depth (Ws(X0)) in [km]: "))
#    dip = float(input("enter the dip angle (dWs/(dWs/dx)(X0)) in [degrees]: "))
#    alpha_ratio = 2 * (Ws / math.tan(math.radians(dip)))
#    return alpha_ratio

# Te estimation from ratio calculation with real depth and has vertical exaggeration
def VE():
    dip = TD.true_dip()
    return
VE()

