import math
def VE_calc():
    H = float(input("enter the length of the vertical scale H in px\n"))
    Y = float(input("enter the real length of the vertical scale Y in km\n"))
    L = float(input("enter the length of the horizontal scale L in px\n")) 
    X = float(input("enter the real length of the Horizontal scale X in km\n"))
    VE = (H/L) / (Y/X)

    print("-------------------------------")
    print("for H = ", H, "Y = ", Y, "for L = ", L,"for X = ", X)
    print("the VE = ", VE )
    print("-------------------------------")
    
VE_calc()
