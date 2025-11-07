
import numpy as np
from cepheid_sim.physics import CepheidParams, phase_to_time_s, radius_Rsun, teff, vrad_kms
from cepheid_sim.photometry import magnitude_V

def main():
    params = CepheidParams(P=10.15, R0=60.0, T0=5800.0, AR=np.array([0.10, 0.02]), AT=np.array([0.03, 0.01]),
                           phi_R=np.array([0.0, 1.2]), phi_T=np.array([0.6, 1.7]), mV_mean=3.9, A_V=0.2)
    phase = np.linspace(0, 1, 500, endpoint=False)
    t = phase_to_time_s(phase, params.P)

    mV = magnitude_V(params, t)
    vr = vrad_kms(params, t)
    R  = radius_Rsun(params, t)
    T  = teff(params, t)

    print(f"mV: mean={mV.mean():.3f}, amp(peak-to-peak)={mV.max()-mV.min():.3f} mag")
    print(f"vr: mean={vr.mean():.3f} km/s, amp={vr.max()-vr.min():.3f} km/s")
    print(f"R : mean={R.mean():.2f} R_sun, amp={R.max()-R.min():.2f} R_sun")
    print(f"T : mean={T.mean():.1f} K,   amp={T.max()-T.min():.1f} K")

if __name__ == '__main__':
    main()
