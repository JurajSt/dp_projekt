import numpy as np
import  warnings
# x = elev uhol (sinus elev uhla)  	(x suradnica)
# SNR = merane data					(y suradnica)
# deg = stupen polynomu
def modelSNRxx(x, SNR, deg):
    poly = np.polyfit(x, SNR, deg)
    model = np.polyval(poly, x)
    residuals = SNR - model
    SNR_vector = residuals - np.mean(residuals)

    return [model, residuals, SNR_vector]

def modelSNR(x, SNR, deg):
    with warnings.catch_warnings():
        warnings.filterwarnings("error")
        try:
            poly = np.polyfit(x, SNR, deg)
            model = np.polyval(poly, x)
            residuals = SNR - model
            SNR_vector = residuals - np.mean(residuals)
            return [model, residuals, SNR_vector]
        except np.RankWarning:
            print "Warning: Polyfit may be poorly conditioned"

    warnings.simplefilter('ignore', np.RankWarning)
    poly = np.polyfit(x, SNR, deg)
    model = np.polyval(poly, x)
    residuals = SNR - model
    SNR_vector = residuals - np.mean(residuals)
    return [model, residuals, SNR_vector]

