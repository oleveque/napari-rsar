from typing import Any, List
from pathlib import Path
import numpy as np

def write_pamela_image(path: str, data: Any, attributes: dict) -> List[str]:
    filename = Path(path)

    # Create .inf file
    infpath = filename.with_suffix('.inf')
    with infpath.open('wt') as inf:
        inf.write(str(data.shape[0]) + '\n')
        inf.write(str(data.shape[1]) + '\n')
        inf.write('1 float\n')
        inf.write('TYPE 1=real\n')
        inf.write("PROGRAMME ORIGINE=napari-rsar plugin")

    # Create .1 file
    slcpath = filename.with_suffix('.1')
    slc = np.memmap(slcpath, mode='w+', dtype=data.dtype, shape=data.shape)
    slc[:] = data[0][:] # Write the first layer of Mipmaps # TODO: Save complex data

    # Return paths to the created files
    return [infpath, slcpath]