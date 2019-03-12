import math
import true_dip_with_vertical_exaggeration as TD

# Te estimation from ratio calculation with real depth
def RD():
    Ws = float(input("enter the depth (Ws(X0)) in [km]: "))
    dip = float(input("enter the dip angle (dWs/(dWs/dx)(X0)) in [degrees]: "))
    alpha_ratio = 2 * (Ws / math.tan(math.radians(dip))) / (1 + math.exp(-math.pi / 2))
    return alpha_ratio, dip

# Te estimation from ratio calculation with real depth and has vertical exaggeration
def VE():
    Ws = float(input("enter the depth (Ws(X0)) in [km]: "))
    dip = TD.true_dip()
    alpha_ratio = 2 * (Ws / math.tan(math.radians(dip)))/ (1 + math.exp(-math.pi / 2))
    return alpha_ratio, dip

# Te estimation from ratio calculation with TWTT assumming piecewise velocity structure:
# from Eldhom et al., 1995 that velocity of sediment is v1 = 3.5 and lava is v2 = 5 [km/s]
def TWTT():
    triangle_x = float(input("enter the triangle x in [km]: "))  # the width for estimate slope
    v_avg1 = 4.;  # V_avg1 is the average velocity for the whole wedge
    v_avg2 = 5.;  #and V_avg2 is the average velocity at Ws(x0) [km/s]
    t_s  = float(input("enter the TWTT(two way travel time) of bottom layer of sediment in sec: ")) 
    t_0  = float(input("enter the time_0 in TWTT in sec: ")) # time_0 is the TWTT for shallower point
    t_1  = float(input("enter the time_1 in TWTT in sec: "))
    ratio =  (v_avg1 / v_avg2) * triangle_x * ((t_0 + t_1)/4. - t_s/2.) / ((t_1 - t_0) / 2.)
    alpha_ratio = ratio * 2/ (1 + math.exp(-math.pi / 2))
    return alpha_ratio
