import numpy as np
from scipy.constants import c as c0
from sargeom.coordinates import Cartographic, CartesianECEF

def find_global_time_offset(data, header):
    # TX position @ integration center
    gp_otx_center = Cartographic(
        longitude=header['data']['log']['flash']['azimuth']['integration']['tx']['center']['position_geo'][0], # [deg]
        latitude=header['data']['log']['flash']['azimuth']['integration']['tx']['center']['position_geo'][1], # [deg]
        height=header['data']['log']['flash']['azimuth']['integration']['tx']['center']['position_geo'][2] # [m]
    ).to_ecef()

    # RX position @ integration center
    gp_orx_center = Cartographic(
        longitude=header['data']['log']['flash']['azimuth']['integration']['rx']['center']['position_geo'][0], # [deg]
        latitude=header['data']['log']['flash']['azimuth']['integration']['rx']['center']['position_geo'][1], # [deg]
        height=header['data']['log']['flash']['azimuth']['integration']['rx']['center']['position_geo'][2] # [m]
    ).to_ecef()

    # Scene center position
    gp_truth = Cartographic(
        longitude=header['data']['log']['flash']['scene']['center_geo'][0], # [deg]
        latitude=header['data']['log']['flash']['scene']['center_geo'][1], # [deg]
        height=header['data']['log']['flash']['scene']['center_geo'][2] # [m]
    ).to_ecef()

    # RX velocity vector @ integration center
    vrx = CartesianECEF(
        x=header['data']['log']['flash']['azimuth']['integration']['rx']['center']['velocity_ecef'][0], # [m/s]
        y=header['data']['log']['flash']['azimuth']['integration']['rx']['center']['velocity_ecef'][0], # [m/s]
        z=header['data']['log']['flash']['azimuth']['integration']['rx']['center']['velocity_ecef'][0] # [m/s]
    )

    # Data axes information
    ax_shape = (header['data']['row']['size'], header['data']['col']['size'])
    ax_origin = (header['data']['row']['origin'], header['data']['col']['origin'])
    ax_step = (header['data']['row']['step'], header['data']['col']['step'])
    ax_name = (header['data']['row']['name'], header['data']['col']['name'])
    ax_unit = (header['data']['row']['unit'], header['data']['col']['unit'])

    # Validate axes information
    assert ax_unit == ('deg', 'deg'), "Units must be degrees."
    assert ax_shape == data.shape, "Data shape must match the shape described in the header."

    # Find indexes of image maximum amplitude
    max_indexes = np.unravel_index(np.argmax(np.abs(data), axis=None), ax_shape)
    if ax_name == ('latitude', 'longitude'):
        longitude_max = ax_origin[1] + max_indexes[1] * ax_step[1]
        latitude_max = ax_origin[0] + max_indexes[0] * ax_step[0]
    elif ax_name == ('longitude', 'latitude'):
        longitude_max = ax_origin[0] + max_indexes[0] * ax_step[0]
        latitude_max = ax_origin[1] + max_indexes[1] * ax_step[1]
    else:
        raise ValueError("Axes names must be 'latitude' and 'longitude'.")
    
    gp_max = Cartographic(
        longitude=longitude_max, # [deg]
        latitude=latitude_max, # [deg]
        height=header['data']['log']['flash']['scene']['center_geo'][-1] # [m]
    ).to_ecef()

    ##############################
    #  Stop-and-Go approximation #
    ##############################
    # delay_truth = (
    #     (gp_otx_center - gp_truth).magnitude() +
    #     (gp_orx_center - gp_truth).magnitude()
    # ) / c0

    # delay_measured = (
    #     (gp_otx_center - gp_max).magnitude() +
    #     (gp_orx_center - gp_max).magnitude()
    # ) / c0

    #####################################
    #  Fast-Time Mean FFL approximation #
    #####################################
    tp_truth = gp_truth - gp_otx_center
    rp_truth = gp_truth - gp_orx_center
    delay_truth = (
        (tp_truth.magnitude() + rp_truth.magnitude()) /
        (c0 + rp_truth.normalize().dot(vrx))
    )

    tp_measured = gp_max - gp_otx_center
    rp_measured = gp_max - gp_orx_center
    delay_measured = (
        (tp_measured.magnitude() + rp_measured.magnitude()) /
        (c0 + rp_measured.normalize().dot(vrx))
    )

    return np.squeeze(delay_measured - delay_truth) # [s]

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

    # Find global time offset
    time_offset = find_global_time_offset(data, header)
    print(time_offset)