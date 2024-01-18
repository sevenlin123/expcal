import numpy as np

readout_time = 20.6
area = 9.7 #m^2
readout = 7 #electrons, readout noise 
pixscale = 0.264
sky_brightness ={ 0:{'u':22.8, 'g':22.1, 'r':21.1, 'i':20.1, 'z':18.7, 'y':18.0},
                  1:{'u':22.6, 'g':22.0, 'r':21.0, 'i':20.1, 'z':18.7, 'y':18.0},
                  2:{'u':22.5, 'g':21.9, 'r':21.0, 'i':20.1, 'z':18.7, 'y':18.0},
                  3:{'u':22.3, 'g':21.8, 'r':20.9, 'i':20.1, 'z':18.7, 'y':18.0},
                  4:{'u':21.9, 'g':21.7, 'r':20.8, 'i':20.1, 'z':18.7, 'y':18.0},
                  5:{'u':21.5, 'g':21.5, 'r':20.8, 'i':20.0, 'z':18.7, 'y':18.0},
                  6:{'u':21.1, 'g':21.4, 'r':20.7, 'i':20.0, 'z':18.6, 'y':18.0},
                  7:{'u':22.7, 'g':21.2, 'r':20.6, 'i':19.9, 'z':18.6, 'y':18.0},
                  8:{'u':20.2, 'g':20.9, 'r':20.5, 'i':19.8, 'z':18.6, 'y':18.0},
                  9:{'u':19.7, 'g':20.7, 'r':20.3, 'i':19.8, 'z':18.5, 'y':17.9},
                 10:{'u':19.2, 'g':20.4, 'r':20.2, 'i':19.7, 'z':18.5, 'y':17.9},
                 11:{'u':18.8, 'g':20.2, 'r':20.1, 'i':19.6, 'z':18.4, 'y':17.9},
                 12:{'u':18.5, 'g':19.9, 'r':20.0, 'i':19.6, 'z':18.4, 'y':17.8},
                 13:{'u':18.1, 'g':19.7, 'r':19.8, 'i':19.5, 'z':18.3, 'y':17.8},
                 14:{'u':17.7, 'g':19.4, 'r':19.7, 'i':19.4, 'z':18.2, 'y':17.7}} #>50 deg from the bright Moon
Flux0 = {'u':7.34e+9, 'g':1.74e+10, 'r':1.23e+10, 'i':1.05e+10, 'z':8.79e+9, 'y':2.76e+9}
throughput = {'u':0.11, 'g':0.33, 'r':0.46, 'i':0.52, 'z':0.52, 'y':0.29} #electrons/photon


def cal_snr_nexp(exptime, total_time=14400, mag = 26.5, airmass = 1, ext_coeff = -0.31, seeing=1, rate=0, Filter='r', moon=7):
    total_nexp = (total_time / (exptime + readout_time))
    total_exp = total_nexp * exptime
    trail_len = rate/3600 * exptime # "/hr  
    mag0 = (total_exp * area * Flux0[Filter]) * 10**(((airmass - 1) * ext_coeff)/2.5) #photon from 0 mag source
    signal0 = mag0*throughput[Filter] #signal from 0 mag source
    signal = signal0*10**(mag/(-2.5)) #signal from specified source
    aperture = np.pi * (2.04*seeing/2)**2 + 2.04*seeing*trail_len
    sky_signal = 10**(sky_brightness[moon][Filter]/-2.5)*aperture*signal0
    ronsq = readout**2*(aperture/(pixscale**2))
    final_snr = signal/(signal+sky_signal+total_nexp*ronsq)**0.5
    return final_snr, total_nexp

def cal_snr(exptime, mag = 26.5, airmass = 1, ext_coeff = -0.31, seeing=1, rate=0, Filter='r', moon=7):
    trail_len = rate/3600 * exptime # "/hr  
    mag0 = (exptime * area * Flux0[Filter]) * 10**(((airmass - 1) * ext_coeff)/2.5) #photon from 0 mag source
    signal0 = mag0*throughput[Filter] #signal from 0 mag source 
    signal = signal0*10**(mag/(-2.5)) #signal from specified source
    aperture = np.pi * (2.04*seeing/2)**2 + 2.04*seeing*trail_len
    sky_signal = 10**(sky_brightness[moon][Filter]/-2.5)*aperture*signal0
    ronsq = readout**2*(aperture/(pixscale**2))
    final_snr = signal/(signal+sky_signal+ronsq)**0.5
    return final_snr

def cal_exp(snr, mag = 26.5, airmass = 1, ext_coeff = -0.31, seeing=1, rate=0, Filter='r', moon=7, niter=10):
    exptime = 1
    for i in range(niter):
        trail_len = rate/3600 * exptime # "/hr  
        mag0 = (exptime * area * Flux0[Filter]) * 10**(((airmass - 1) * ext_coeff)/2.5) #photon from 0 mag source
        signal0 = mag0*throughput[Filter] #signal from 0 mag source 
        aperture = np.pi * (2.04*seeing/2)**2 + 2.04*seeing*trail_len
        sky_signal = 10**(sky_brightness[moon][Filter]/-2.5)*aperture*signal0
        ronsq = readout**2*(aperture/(pixscale**2))
        final_s = (snr**2 + snr*(snr**2+4*(sky_signal+ronsq))**0.5)/2
        final_s0 = final_s * 10**(mag/2.5)
        exptime*= final_s0/signal0

    return exptime