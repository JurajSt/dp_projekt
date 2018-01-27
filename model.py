import numpy as np

# x = elev uhol (sinus elev uhla)  	(x suradnica)
# SNR = merane data					(y suradnica)
# deg = stupen polynomu
def modelSNR(x, SNR, deg):
    poly = np.polyfit(x, SNR, deg)
    model = np.polyval(poly, x)
    residuals = SNR - model
    #SNR_vector = residuals - np.mean(residuals)

    return [model, residuals] #SNR_vector]