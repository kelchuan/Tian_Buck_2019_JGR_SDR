import math
import true_dip_with_vertical_exaggeration as TD
import Df_Te_alpha as DfTe
def main():
    Depth = float(input("enter the depth (Ws(X0)) in [km]: "))
    dip = TD.true_dip()
    alpha = 2 * (Depth / math.tan(math.radians(dip)))
    alpha = alpha * 1000
    
    kappa,alpha_ratio,Te_ratio = DfTe.Df_Te()
    Te = (alpha / kappa)**(4.0/3.0) #[m]
    print ("Depth is:", Depth, "[m]", "dip angle is:", dip, "[m]", "and Te_ratio is:", Te, "[m]", "alpha_ratio is ", alpha, "[m]")



main()
