import math
import true_dip_with_vertical_exaggeration as TD
import Df_Te_alpha as DfTe
def main():
    Depth = float(input("enter the depth (Ws(X0)) in [km]: "))
    dip = TD.true_dip()
    alpha = 2 * (Depth / math.tan(math.radians(dip)))
    alpha = alpha * 1000
    
    kappa = DfTe.Df_Te()
    Te = (alpha / kappa)**(4.0/3.0) #[m]
    print ("Depth is:", Depth, "[m]", "dip angle is:", dip, "[m]", "and Te_ratio is:", Te, "[m]", "alpha_ratio is ", alpha, "[m]")
    print ("__________________")
    print ("Results from ratio estimation")
    triangle_x = float(input("enter the triangle x in [km]: ")) 
    v1 = 3.5; v2 = 5.; # V1 is the sediment velocity and V2 is the lava velocity 
    t_s  = float(input("enter the sediment thickness in TWTT in sec: ")) 
    t_o  = float(input("enter the T_o in TWTT in sec: "))
    t_1  = float(input("enter the T_1 in TWTT in sec: ")) 
    
    alpha_ratio = (ts * 2. / (t_1 - t_o)) * (v1 / v2) * triangle_x + ((t_o + t_1 - 2*t_s)* triangle_x / 2.) / ((t_1 - t_o) / 2)
    
    Te_ratio = (alpha_ratio / kappa)**(4.0/3.0) #[m]
    print (alpha_ratio, Te_ratio)  
main()
