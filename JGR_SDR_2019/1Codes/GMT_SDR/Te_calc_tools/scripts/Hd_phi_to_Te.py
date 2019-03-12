import math
import Xf_Te as XfTe
import Ratio_Te as RTe
def main():
    kappa,alpha_xf,Te_xf,coeff = XfTe.Xf_Te()
    hd = float(input("input hd in meters\n"))
    phi = float(input("input phi in degrees\n"))
    print(coeff,hd,phi)
    Te = (hd / math.tan(math.radians(phi)) * coeff)**(4./3.)
    print("Te_phi=",Te)
main()
