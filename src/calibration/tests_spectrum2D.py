
import numpy as np
from scipy.constants import c as c0
from sargeom.coordinates import Cartographic, CartesianECEF, Cartesian3

from matplotlib import pyplot as plt

ZAXIS = Cartesian3.UNIT_Z()

# def ground_bisector_vector_SG(otx, orx):
#     txp,  rxp = -otx, -orx
#     utxp, urxp = txp.normalize(), rxp.normalize()
#     beta = utxp + urxp # Bisector vector using Stop-And-Go approximation
#     return beta.reject_from(ZAXIS) # Ground bisector vector

def ground_bisector_vector_FFL(otx, orx, vrx):
    txp,  rxp = -otx, -orx
    utxp, urxp = txp.normalize(), rxp.normalize()
    beta = c0 * (utxp + urxp) / (c0 + urxp.dot(vrx)) # Bisector vector using the Far-Field Linear approximation
    return beta.reject_from(ZAXIS) # Ground bisector vector

def spectrum2(data, header):
    # Scene center position
    gp_center = Cartographic(
        longitude=header['data']['log']['flash']['scene']['center_geo'][0], # [deg]
        latitude=header['data']['log']['flash']['scene']['center_geo'][1], # [deg]
        height=header['data']['log']['flash']['scene']['center_geo'][2] # [m]
    )

    # TX integration center position
    otx_center = Cartographic(
        longitude=header['data']['log']['flash']['azimuth']['integration']['tx']['center']['position_geo'][0], # [deg]
        latitude=header['data']['log']['flash']['azimuth']['integration']['tx']['center']['position_geo'][1], # [deg]
        height=header['data']['log']['flash']['azimuth']['integration']['tx']['center']['position_geo'][2] # [m]
    ).to_ecef().to_enu(origin=gp_center)

    # RX integration center position
    orx_center = Cartographic(
        longitude=header['data']['log']['flash']['azimuth']['integration']['rx']['center']['position_geo'][0], # [deg]
        latitude=header['data']['log']['flash']['azimuth']['integration']['rx']['center']['position_geo'][1], # [deg]
        height=header['data']['log']['flash']['azimuth']['integration']['rx']['center']['position_geo'][2] # [m]
    ).to_ecef().to_enu(origin=gp_center)

    # Velocity vectors estimation @ integration center
    vrx_center = CartesianECEF(
        x=header['data']['log']['flash']['azimuth']['integration']['rx']['center']['velocity_ecef'][0], # [m/s]
        y=header['data']['log']['flash']['azimuth']['integration']['rx']['center']['velocity_ecef'][1], # [m/s]
        z=header['data']['log']['flash']['azimuth']['integration']['rx']['center']['velocity_ecef'][2] # [m/s]
    ).to_enuv(origin=gp_center)
    
    # Ground bisector vector @ scene center
    betag_center = ground_bisector_vector_FFL(otx_center, orx_center, vrx_center)

    # Wavenumber @ integration center min and max (relative to frequency)
    fc = header['data']['log']['range']['center_frequency_hz']
    bproc = header['data']['log']['range']['processed_bandwidth_hz']
    kc_eta_center   = 2 * np.pi / c0 * fc * betag_center
    kmin_eta_center = 2 * np.pi / c0 * (fc - 0.5 * bproc) * betag_center - kc_eta_center
    kmax_eta_center = 2 * np.pi / c0 * (fc + 0.5 * bproc) * betag_center - kc_eta_center

    # TX integration start position
    otx_start = Cartographic(
        longitude=header['data']['log']['flash']['azimuth']['integration']['tx']['start']['position_geo'][0], # [deg]
        latitude=header['data']['log']['flash']['azimuth']['integration']['tx']['start']['position_geo'][1], # [deg]
        height=header['data']['log']['flash']['azimuth']['integration']['tx']['start']['position_geo'][2] # [m]
    ).to_ecef().to_enu(origin=gp_center)

    # RX integration center position
    orx_start = Cartographic(
        longitude=header['data']['log']['flash']['azimuth']['integration']['rx']['start']['position_geo'][0], # [deg]
        latitude=header['data']['log']['flash']['azimuth']['integration']['rx']['start']['position_geo'][1], # [deg]
        height=header['data']['log']['flash']['azimuth']['integration']['rx']['start']['position_geo'][2] # [m]
    ).to_ecef().to_enu(origin=gp_center)

    # Velocity vectors estimation @ integration start
    vrx_start = CartesianECEF(
        x=header['data']['log']['flash']['azimuth']['integration']['rx']['start']['velocity_ecef'][0], # [m/s]
        y=header['data']['log']['flash']['azimuth']['integration']['rx']['start']['velocity_ecef'][1], # [m/s]
        z=header['data']['log']['flash']['azimuth']['integration']['rx']['start']['velocity_ecef'][2] # [m/s]
    ).to_enuv(origin=gp_center)
    
    # Ground bisector vector @ start
    betag_start = ground_bisector_vector_FFL(otx_start, orx_start, vrx_start)
    
    # Wavenumber @ integration start min and max (relative to frequency)
    kmin_eta_start = 2 * np.pi / c0 * (fc - 0.5 * bproc) * betag_start - kc_eta_center
    kmax_eta_start = 2 * np.pi / c0 * (fc + 0.5 * bproc) * betag_start - kc_eta_center

    # TX integration end position
    otx_end = Cartographic(
        longitude=header['data']['log']['flash']['azimuth']['integration']['tx']['end']['position_geo'][0], # [deg]
        latitude=header['data']['log']['flash']['azimuth']['integration']['tx']['end']['position_geo'][1], # [deg]
        height=header['data']['log']['flash']['azimuth']['integration']['tx']['end']['position_geo'][2] # [m]
    ).to_ecef().to_enu(origin=gp_center)

    # RX integration center position
    orx_end = Cartographic(
        longitude=header['data']['log']['flash']['azimuth']['integration']['rx']['end']['position_geo'][0], # [deg]
        latitude=header['data']['log']['flash']['azimuth']['integration']['rx']['end']['position_geo'][1], # [deg]
        height=header['data']['log']['flash']['azimuth']['integration']['rx']['end']['position_geo'][2] # [m]
    ).to_ecef().to_enu(origin=gp_center)

    # Velocity vectors estimation @ integration end
    vrx_end = CartesianECEF(
        x=header['data']['log']['flash']['azimuth']['integration']['rx']['end']['velocity_ecef'][0], # [m/s]
        y=header['data']['log']['flash']['azimuth']['integration']['rx']['end']['velocity_ecef'][1], # [m/s]
        z=header['data']['log']['flash']['azimuth']['integration']['rx']['end']['velocity_ecef'][2] # [m/s]
    ).to_enuv(origin=gp_center)
    
    # Ground bisector vector @ end
    betag_end = ground_bisector_vector_FFL(otx_end, orx_end, vrx_end)
    
    # Wavenumber @ integration end min and max (relative to frequency)
    kmin_eta_end = 2 * np.pi / c0 * (fc - 0.5 * bproc) * betag_end - kc_eta_center
    kmax_eta_end = 2 * np.pi / c0 * (fc + 0.5 * bproc) * betag_end - kc_eta_center

    # Image spacing @ scene center
    spacing_x_m = header['data']['log']['flash']['scene']['spacing_x_m']
    spacing_y_m = header['data']['log']['flash']['scene']['spacing_y_m']

    # Wavenumber extent
    kxmin, kxmax = -np.pi / spacing_x_m, np.pi / spacing_x_m
    kymin, kymax = -np.pi / spacing_y_m, np.pi / spacing_y_m

    fft_data = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(data)))

    # see: https://stackoverflow.com/a/56595416

    black = np.full((*fft_data.shape, 4), 0.)
    black[:,:,-1] = np.abs(fft_data) / np.abs(fft_data).max()
    black[:,:,-1] = 1 - black[:,:,-1]

    plt.figure()
    plt.imshow(
        np.angle(fft_data), cmap='hsv',
        extent=[kxmin, kxmax, kymin, kymax],
        interpolation='nearest'
    )
    plt.imshow(
        black,
        extent=[kxmin, kxmax, kymin, kymax],
        interpolation='nearest'
    )
    plt.plot(
        [kmin_eta_start.x, kmax_eta_start.x],
        [kmin_eta_start.y, kmax_eta_start.y],
        '-o', color='tab:red',
        label='@Integration start'
    )
    plt.plot(
        [kmin_eta_center.x, kmax_eta_center.x],
        [kmin_eta_center.y, kmax_eta_center.y],
        '-o', color='tab:green',
        label='@Integration center'
    )
    plt.plot(
        [kmin_eta_end.x, kmax_eta_end.x],
        [kmin_eta_end.y, kmax_eta_end.y],
        '-o', color='tab:blue',
        label='@Integration end'
    )
    plt.xlabel("kx ground wavenumber [rad/m]")
    plt.ylabel("ky ground wavenumber [rad/m]")
    plt.legend()
    plt.grid()
    plt.show()

if __name__ == '__main__':
    from pathlib import Path
    import tomlkit
    import sys

    # Get filename
    filename = Path(sys.argv[1])

    # Load data and header
    with open(filename, 'r') as f:
        header = tomlkit.load(f)
        data = np.memmap(filename.with_suffix(''),
            shape=(header['data']['row']['size'], header['data']['col']['size']),
            dtype='c8'
        )

    # Plot spectrum
    spectrum2(data, header)