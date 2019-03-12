import math
def true_dip():
    VE = float(input("enter the Vertical Exaggeration: (Enter 0 if you have scale rather than VE available)"))

    if VE == 0.0:
        horizontal = float(input("enter horizontal scale: "))
        vertical = float(input("enter corresponding vertical scale(same apparent length of horizontal scale): "))
        VE = horizontal/vertical

    dip_VE = float(input("enter the apparent dip angle: "))
    dip = math.degrees(math.atan(1/VE * math.tan(math.radians(dip_VE))))
    
    print ("The true dip angle is:", dip, "with apparent dip of:", dip_VE, "and VE of ", VE)
    return dip

