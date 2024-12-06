import numpy as np
from scipy.constants import c as c0
from sargeom.coordinates import Cartographic, CartesianECEF

def find_global_time_offset(data, header):
    # RX integration center position
    gp_orx_center = Cartographic(
        longitude=header['data']['log']['flash']['azimuth']['integration']['rx']['center']['position_geo'][0], # [deg]
        latitude=header['data']['log']['flash']['azimuth']['integration']['rx']['center']['position_geo'][1], # [deg]
        height=header['data']['log']['flash']['azimuth']['integration']['rx']['center']['position_geo'][2] # [m]
    ).to_ecef()

    # TX integration center position
    gp_otx_center = Cartographic(
        longitude=header['data']['log']['flash']['azimuth']['integration']['tx']['center']['position_geo'][0], # [deg]
        latitude=header['data']['log']['flash']['azimuth']['integration']['tx']['center']['position_geo'][1], # [deg]
        height=header['data']['log']['flash']['azimuth']['integration']['tx']['center']['position_geo'][2] # [m]
    ).to_ecef()

    # Ground truth position
    gp_truth = Cartographic(
        longitude=header['data']['log']['flash']['scene']['center_geo'][0], # [deg]
        latitude=header['data']['log']['flash']['scene']['center_geo'][1], # [deg]
        height=header['data']['log']['flash']['scene']['center_geo'][2] # [m]
    ).to_ecef()

    # Velocity vectors estimation @ integration center
    vrx = CartesianECEF(
        x=header['data']['log']['flash']['azimuth']['integration']['rx']['center']['velocity_ecef'][0], # [m/s]
        y=header['data']['log']['flash']['azimuth']['integration']['rx']['center']['velocity_ecef'][0], # [m/s]
        z=header['data']['log']['flash']['azimuth']['integration']['rx']['center']['velocity_ecef'][0] # [m/s]
    )
    vtx = CartesianECEF(
        x=header['data']['log']['flash']['azimuth']['integration']['tx']['center']['velocity_ecef'][0], # [m/s]
        y=header['data']['log']['flash']['azimuth']['integration']['tx']['center']['velocity_ecef'][1], # [m/s]
        z=header['data']['log']['flash']['azimuth']['integration']['tx']['center']['velocity_ecef'][2] # [m/s]
    )

    # 
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
    #     (gp_orx_center - gp_truth).norm() +
    #     (gp_otx_center - gp_truth).norm()
    # ) / c0

    # delay_measured = (
    #     (gp_orx_center - gp_max).norm() +
    #     (gp_otx_center - gp_max).norm()
    # ) / c0

    #####################################
    #  Fast-Time Mean FFL approximation #
    #####################################
    rp_truth = gp_truth - gp_orx_center
    tp_truth = gp_truth - gp_otx_center
    delay_truth = (
        (rp_truth.norm() + tp_truth.norm()) /
        (c0 + rp_truth.normalize().dot(vrx))
    )

    rp_measured = gp_max - gp_orx_center
    tp_measured = gp_max - gp_otx_center
    delay_measured = (
        (rp_measured.norm() + tp_measured.norm()) /
        (c0 + rp_measured.normalize().dot(vrx))
    )

    return delay_measured - delay_truth # [s]