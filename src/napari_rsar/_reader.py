from tomlkit import load as parse
from pathlib import Path
import numpy as np
import napari

def rsar_file_reader(path):
    """
    Reads rSAR header information and data from the specified file.

    Parameters
    ----------
    path : str
        Path to the rSAR header file.

    Returns
    -------
    list
        List containing LayerData tuples.
    """
    filename = Path(path)

    # Read header information from the specified file
    with filename.open('rt', encoding='utf-8') as h:
        header = parse(h)

    # Check if the header information is valid
    assert header['data']['version'] == 1, "Invalid version number in the header information"
    assert header['data']['kind'] in ['slc', 'rcc', 'dem', 'lobe'], "Invalid data kind in the header information"

    # Determine the byte order based on the header information
    byteorder = '<' if header['data']['endianness'] == 'le' else '>'

    # Determine the data type based on the header information
    datatype = header['data']['datatype']
    if datatype == 'u16':
        if header['data']['is_complex']:
            raise ValueError("Unsupported data type")
        else:
            type = 'f2' # 16-bit floating-point number
    elif datatype == 'f32':
        if header['data']['is_complex']:
            type = 'c8' # 64-bit complex floating-point number
        else:
            type = 'f4' # 34-bit floating-point number
    else:
        raise ValueError("Unsupported data type")

    # Load the data from the file using numpy memmap
    data = np.memmap(
        filename.with_suffix(''),
        shape=(header['data']['row']['size'], header['data']['col']['size']),
        offset=header['data']['global_offset_byte'],
        dtype=byteorder+type,
        mode='r'
    )

    # Extract scale and offset values for data rescaling
    scale = header['data']['data']['scale']
    offset = header['data']['data']['offset']
    
    if header['data']['is_complex']:
        # Generate mipmaps for the image data
        index = 0
        layer = scale * (data + offset)
        layer1_mipmap = [np.abs(layer)]
        layer2_mipmap = [np.angle(layer)]
        while max(layer1_mipmap[-1].shape) > 1024:
            new_data = layer1_mipmap[index][::2, ::2]
            layer1_mipmap.append(new_data)

            new_data = layer2_mipmap[index][::2, ::2]
            layer2_mipmap.append(new_data)

            index += 1

        # Create a layer for the magnitude data
        (data1, attributes1, type1) = napari.layers.Image(
            data=layer1_mipmap,
            name=f'[abs] {filename.stem}',
            contrast_limits=[0, np.percentile(layer1_mipmap[0], 92)],
            axis_labels=(header['data']['row']['name'], header['data']['col']['name']),
            interpolation2d='nearest',
            multiscale=True,
            colormap='gray',
            visible=True,
            metadata=header['data']
        ).as_layer_data_tuple()

        # Create a second layer for the angle data
        (data2, attributes2, type2) = napari.layers.Image(
            data=layer2_mipmap,
            name=f'[ang] {filename.stem}',
            contrast_limits = [-np.pi, np.pi],
            axis_labels=(header['data']['row']['name'], header['data']['col']['name']),
            interpolation2d='nearest',
            multiscale=True,
            colormap='hsv',
            visible=False,
            metadata=header['data']
        ).as_layer_data_tuple()
        
        return [(data2, dict(attributes2), type2), (data1, dict(attributes1), type1)]
    
    else:
        # Generate mipmaps for the image data
        index = 0
        layer_mipmap = [scale * (data + offset)]
        while max(layer_mipmap[-1].shape) > 1024:
            new_data = layer_mipmap[index][::2, ::2]
            layer_mipmap.append(new_data)
            index += 1
        
        # Create a layer for the image data
        (data, attributes, type) = napari.layers.Image(
            data=layer_mipmap,
            name=filename.stem,
            contrast_limits=[0, np.percentile(layer_mipmap[0], 92)],
            axis_labels=(header['data']['row']['name'], header['data']['col']['name']),
            interpolation2d='nearest',
            multiscale=True,
            colormap='gray',
            visible=True,
            metadata=header['data']
        ).as_layer_data_tuple()
        
        return [(data, dict(attributes), type)]

def get_rsar_reader(path):
    """
    Returns a rSAR file reader function if the given path is a valid rSAR file.

    Parameters
    ----------
    path : str
        Path to the rSAR file.

    Returns
    -------
    callable or None
        If the path is valid, returns the rSAR file reader function. Otherwise, returns None.
    """
    suffixes = ['.slc.toml', '.rcc.toml', '.dem.toml', '.lobe.toml']

    # Check if the path is a string and ends with any of the specified suffixes
    if isinstance(path, str) and any(path.endswith(ext) for ext in suffixes):
        return rsar_file_reader # Return the rSAR file reader function
    else:
        return None # Return None if the path is not valid