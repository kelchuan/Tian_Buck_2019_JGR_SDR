import math
import Xf_Te as XfTe
import Ratio_Te as RTe
# This main function will do two things:
# 1. call Xf_Te function to init parameters and to get kappa, alpha_xf and Te_xf
# 2. call ratio_Te function to get alpha_ratio and Te_ratio
  # 2.1 determine whether data is twtt or real depth or real depth with V.E.
def main():
    # The main function has two parts, one is Te estimation from Xf, the other is Te estimation from ratio
    # 1. Te estimation from Xf
    kappa,alpha_xf,Te_xf,coeff = XfTe.Xf_Te()
    # 2. Te estimation from ratio of Ws/(dWs/dx)
    data_type = input("Twtt(tw) or real depth(rd) or real depth with vertical exageration(ve)? (tw/rd/ve)");
    if data_type == 'tw':
        alpha_ratio = RTe.TWTT()
    elif data_type == 'rd':
        alpha_ratio,phi = RTe.RD()
    elif data_type == 've':
        alpha_ratio,phi = RTe.VE()
    alpha_ratio = alpha_ratio * 1000 #[m]
    Te_ratio = (alpha_ratio / kappa)**(4.0/3.0) #[m]

    print ( "Te_xf is:", Te_xf, "[m]", "alpha_xf is ", alpha_xf, "[m]")
    print ( "Te_ratio is:", Te_ratio, "[m]", "alpha_ratio is ", alpha_ratio, "[m]")
    #    print ("Depth is:", Depth, "[m]", "dip angle is:", dip, "[m]", "and Te_ratio is:", Te, "[m]", "alpha_ratio is ", alpha, "[m]")
    # error estimation
#    err_abs = math.fabs(Te_xf - Te_ratio)
    err_abs = Te_xf - Te_ratio
    err_percent = err_abs / Te_xf * 100
    print("-----------------------------------------------------")
    print("For Te estimation, absolute error is ", err_abs, "relative error compared to result from Xf is \n Te_err_percentage = ", err_percent, "%")
    print("-----------------------------------------------------")

#    err_abs_alpha = math.fabs(alpha_xf - alpha_ratio)
    err_abs_alpha = alpha_xf - alpha_ratio
    err_percent_alpha = err_abs_alpha / alpha_xf * 100
    print("-----------------------------------------------------")
    print("For alpha estimation, absolute error is ", err_abs_alpha, "relative error compared to result from Xf is \n Alpha_err_percentage", err_percent_alpha, "%")
    print("-----------------------------------------------------")

    alpha_avg = 0.5 * (alpha_xf + alpha_ratio)
    Te_avg = 0.5 * (Te_xf + Te_ratio)

    hd = math.tan(math.radians(phi)) / coeff * Te_avg**0.75
    print("-----------------------------------------------------")
    print("The estimated hd from phi and calculated Te_xf, Te_gamma is\n Hd = %.2f m" % hd)
    print("-----------------------------------------------------")
    if err_percent_alpha < 20:
        print("-----------------------------------------------------")
        print("The estimations are within 20% accuracy!!")
        print("The relative error in alpha is %.2f" % err_percent_alpha, "%")
        print("The relative error in Te is (from Xf compared to ratio)\n Te_err_percentage = %.2f" % err_percent, "%")
#        print("The average Te is", Te_avg)
        print("The average Te is %.2f" % Te_avg)
#        print("The average alpha is ",  alpha_avg)
        print("The average alpha is %.2f" % alpha_avg)
        print("-----------------------------------------------------")
    else:
        print("-----------------------------------------------------")
        print("The estimations are NOT within 20% accuracy!!")
        print("The relative error in Te is (from Xf compared to ratio)\n Te_err_percentage = %.2f" % err_percent, "%")
#        print("The average Te is", Te_avg)
        print("The average Te is %.2f" % Te_avg)
#        print("The average alpha is ",  alpha_avg)
        print("The average alpha is %.2f" % alpha_avg)       
        print("-----------------------------------------------------")

    
main()
