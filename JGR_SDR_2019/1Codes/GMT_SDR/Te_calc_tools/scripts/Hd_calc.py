import math
import Xf_Te as XfTe
import Ratio_Te as RTe
def main():
    kappa,alpha_xf,Te_xf,coeff = XfTe.Xf_Te()
    Te_avg = float(input("input Te in meters\n"))
    alpha_ratio,phi = RTe.RD()
    hd = math.tan(math.radians(phi)) / coeff * Te_avg**0.75
    print("hd=",hd)
main()
