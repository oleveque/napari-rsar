import numpy as np
from scipy.constants import c as c0
from scipy.interpolate import RegularGridInterpolator
from sargeom.coordinates import Cartographic, CartesianECEF, Cartesian3

from pathlib import Path
from pandas import DataFrame
from matplotlib import pyplot as plt

ZAXIS = Cartesian3.UNIT_Z

HEADER = (
"# Fields description:\n"
"# -------------------\n"
"#    o Frequency axis\n"
"#        - FREQUENCY_HZ [Hz]: The RF frequency range of the estimated waveform correction file in Hertz.\n"
"#    o Complex waveform correction values:\n"
"#        - AMPLITUDE_LIN [-]: The normalized amplitude of the complex waveform correction file.\n"
"#        - PHASE_RAD [rad]: The unwrapped phase of the complex waveform correction file in radians.\n\n")

def estimate_waveform_correction(data, header, output_path, display=False):
    # RX integration center position
    gp_orx_center = Cartographic(
        longitude=header['data']['log']['flash']['azimuth']['integration']['rx']['center']['position_geo'][0], # [deg]
        latitude=header['data']['log']['flash']['azimuth']['integration']['rx']['center']['position_geo'][1], # [deg]
        height=header['data']['log']['flash']['azimuth']['integration']['rx']['center']['position_geo'][2] # [m]
    )

    # TX integration center position
    gp_otx_center = Cartographic(
        longitude=header['data']['log']['flash']['azimuth']['integration']['tx']['center']['position_geo'][0], # [deg]
        latitude=header['data']['log']['flash']['azimuth']['integration']['tx']['center']['position_geo'][1], # [deg]
        height=header['data']['log']['flash']['azimuth']['integration']['tx']['center']['position_geo'][2] # [m]
    )

    # Velocity vectors estimation @ integration center
    vrx = CartesianECEF(
        x=header['data']['log']['flash']['azimuth']['integration']['rx']['center']['velocity_ecef'][0], # [m/s]
        y=header['data']['log']['flash']['azimuth']['integration']['rx']['center']['velocity_ecef'][1], # [m/s]
        z=header['data']['log']['flash']['azimuth']['integration']['rx']['center']['velocity_ecef'][2] # [m/s]
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
    )
    
    # 
    orx_center = gp_orx_center.to_ecef().to_enu(origin=gp_max)
    otx_center = gp_otx_center.to_ecef().to_enu(origin=gp_max)
    vrx = vrx.to_enu(origin=gp_max)
    vtx = vtx.to_enu(origin=gp_max)

    # Bisector vector @ maximum
    txp,  rxp = -otx_center, -orx_center
    utxp, urxp = txp.normalize(), rxp.normalize()
    beta =  c0 * (utxp + urxp) / (c0 + urxp.dot(vrx))

    # Ground bisector vector @ scene center
    betag = beta.reject_from(ZAXIS)

    # First derivative of betagc
    dbetarx = -(vrx - (urxp.dot(vrx)) * urxp) / rxp.norm()
    dbetatx = -(vtx - (utxp.dot(vtx)) * utxp) / txp.norm()
    dbeta = (
        (c0 * (dbetatx + dbetarx) - (dbetarx.dot(vrx)) * beta) /
        (c0 + urxp.dot(vrx))
    )

    # Ground bisector vector @ scene center
    dbetag = dbeta.reject_from(ZAXIS)

    # Ground range direction @ scene center
    rg = dbetag.cross(ZAXIS).normalize()
    rg *= np.sign(rg.dot(betag))
        
    # Ground doppler direction @ scene center
    dg = betag.cross(ZAXIS).normalize()
    dg *= -np.sign(dg.dot(dbetag))
    # note: la direction de "dbeta" est opposée à celle de "v" d'où une direction
    #       d'accroissement des Doppler "inverse" par rapport au sens de dbeta.
    
    # Image spacing @ scene center
    spacing_x_m = header['data']['log']['flash']['scene']['spacing_x_m']
    spacing_y_m = header['data']['log']['flash']['scene']['spacing_y_m']
    
    # Ground extent centered around max_indexes
    xaxis = (np.arange(ax_shape[1]) - max_indexes[1]) * spacing_x_m
    yaxis = (np.arange(ax_shape[0]) - max_indexes[0]) * -spacing_y_m # note: top start at max value
    
    # Create a 2D linear interpolator of the spectrum
    interpolator = RegularGridInterpolator(
        (yaxis, xaxis), data,
        method='linear'
    ) # note: `RegularGridInterpolator` uses 'ij' indexing

    # Interpolation conditions
    dist_max = np.min( # distance in the range/doppler direciton for interpolation
        (ax_shape[1]//2 - np.abs(ax_shape[1]//2 - max_indexes[1])) * 2.0 * spacing_x_m,
        (ax_shape[0]//2 - np.abs(ax_shape[0]//2 - max_indexes[0])) * 2.0 * spacing_y_m,
    )
    spacing_min = np.min(spacing_x_m, spacing_y_m) # taille pixel minimum
    ninterp = int(np.round(dist_max / spacing_min)) + 1 # nombre de points pour l'interpolation
    if ninterp % 2 == 0: # Assure d'avoir un nombre impair (et donc un vrai zéro au centre)
        ninterp -= 1

    # Axe "distance algébrique" d'interpolation
    range_dist = (np.arange(ninterp) - ninterp // 2) * spacing_min
    range_axis = range_dist * rg.expand_dims(1) # Vecteur d'interpolation

    # Interpolation selon l'axe voulu
    slc_range = interpolator((range_axis.y, range_axis.x))
    slc_range /= np.abs(slc_range).max()

    # Range spectrum
    spec_range = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(slc_range)))

    # krmin, krmax = -np.pi / spacing_min, np.pi / spacing_min
    dk = 2 * np.pi / ((ninterp - 1) * spacing_min)
    krange = (np.arange(ninterp) - ninterp // 2) * dk # Nombre d'onde (centré) du spectre distance interpolé

    # Range spectrum interpolation
    fc = header['data']['log']['range']['center_frequency_hz']
    bproc = header['data']['log']['range']['processed_bandwidth_hz']
    df = c0 / (2 * np.pi) * dk / betag.norm() # échantillonnage fréquentiel correspondant à l'échantillonnage en nombre d'onde
    nf = int(bproc / df) + 1
    if nf < 101: # On prend au moins 101 points d'interpolation (notamment si l'image d'origine n'a pas été suréchantillonnée au calcu... Pas bien !!)
        nf = 101
    # df_interp = bproc / (nf - 1)
    finterp = np.linspace(fc-0.5*bproc, fc+0.5*bproc, nf) # Fréquences (RF) d'interpolation
    kinterp = 2 * np.pi / c0 * (finterp - fc) * betag.norm() # Nombre d'onde (centré) correspondant aux fréquences
    
    # Interpolated range spectrum
    spec_range_interp = np.interp(kinterp, krange, spec_range)
    spec_range_interp /= np.abs(spec_range_interp).mean()

    ############################
    # Waveform correction file #
    ############################
    band_path = Path(output_path).with_suffix('.band.csv')
    band_data = DataFrame({
        'FREQUENCY_HZ': finterp,
        'AMPLITUDE_LIN': np.abs(spec_range_interp),
        'PHASE_RAD': np.unwrap(np.angle(spec_range_interp))
    })

    with band_path.open('wt', encoding='utf-8') as f:
        f.write(f"# File: `{band_path.name}`\n\n{HEADER}")
        band_data.to_csv(f, index=False, sep=';', encoding='utf-8', lineterminator='\n')

    ###################
    # Display results #
    ###################
    if display:
        # Axe "azimut algébrique" d'interpolation
        azi_dist = (np.arange(ninterp) - ninterp // 2) * spacing_min
        azi_axis = azi_dist * dg.expand_dims(1) # Vecteur d'interpolation
        
        # Interpolation selon l'axe voulu
        slc_azi = interpolator((azi_axis.y, azi_axis.x))
        slc_azi /= np.abs(slc_azi).max()
        
        #########
        # PLOTS #
        #########
        amp = np.abs(data)
        amp /= np.max(amp)
        xaxis_plot = (np.arange(ax_shape[1]) - ax_shape[1]//2) * spacing_x_m
        yaxis_plot = (np.arange(ax_shape[0]) - ax_shape[0]//2) * -spacing_y_m # note: top start at max value
        center = Cartesian3(
            x=(max_indexes[1] - ax_shape[1]//2) * spacing_x_m, # [m]
            y=(max_indexes[0] - ax_shape[0]//2) * -spacing_y_m, # [m]
            z=0 # [m]
        )
        range_axis_plot = range_axis + center
        azi_axis_plot = azi_axis + center
        fig = plt.figure()
        plt.imshow(
            amp, cmap='gray', interpolation='bilinear', vmax=0.5,
            extent=[xaxis_plot.min(), xaxis_plot.max(),
                    yaxis_plot.min(), yaxis_plot.max()],
        )
        plt.arrow(
            range_axis_plot[0].x, range_axis_plot[0].y,
            range_axis_plot[-1].x - range_axis_plot[0].x,
            range_axis_plot[-1].y - range_axis_plot[0].y,
            length_includes_head=True,
            head_width=dist_max*0.025,
            color='tab:red',
            # label="Ground range"
        )
        plt.arrow(
            azi_axis_plot[0].x, azi_axis_plot[0].y,
            azi_axis_plot[-1].x - azi_axis_plot[0].x,
            azi_axis_plot[-1].y - azi_axis_plot[0].y,
            length_includes_head=True,
            head_width=dist_max*0.025,
            color='tab:blue',
            # label="Ground azimuth"
        )
        plt.plot(
            center.x,
            center.y,
            'o', color='tab:orange'
        )
        plt.xlabel('East ground range [m] (relative to center)')
        plt.ylabel('North ground range [m] (relative to center)')
        plt.axis('equal')
        plt.grid(True)
        plt.tight_layout()
        
        fig, (ax1, ax2) = plt.subplots(2)

        ax1.plot(range_dist, 20*np.log10(np.abs(slc_range)), color='tab:red')
        ax1.set_title('Ground range slice')
        ax1.set_ylabel("Intensity [dB]")
        ax1.set_xlabel("Ground range [m] (relative to center)")
        ax1.grid(True)

        ax2.plot(range_dist, 20*np.log10(np.abs(slc_azi)), color='tab:blue')
        ax2.set_title('Ground azimuth slice')
        ax2.set_xlabel("Ground azimuth [m] (relative to center)")
        ax2.set_ylabel("Intensity [dB]")
        ax2.sharex(ax1)
        ax2.grid(True)
        plt.tight_layout()
        
        fig, (ax1, ax2) = plt.subplots(2)
        fig.suptitle('Range spectrum slice')

        ax1.plot(finterp*1e-9, 10*np.log10(np.abs(spec_range_interp)))
        ax1.set_ylabel("Amplitude [dB]")
        ax1.grid(True)

        ax2.plot(finterp*1e-9, np.rad2deg(np.unwrap(np.angle(spec_range_interp))))
        ax2.set_xlabel("Frequency [GHz]")
        ax2.set_ylabel("Phase [°]")
        ax2.sharex(ax1)
        ax2.grid(True)
        plt.tight_layout()

        # from pymela.visualisation import pqi
        # pqi(slc_range, title='Ground range slice', oversample=1, dist_step=spacing_min)
        # pqi(slc_azi, title='Ground azimuth slice', oversample=1, dist_step=spacing_min)
        # plt.show()