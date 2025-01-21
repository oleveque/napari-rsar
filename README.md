# rsar

A repository that provides a Python API for the rSAR format and a plugin for napari to facilitate the visualization and analysis of rSAR files.

## Installation

1. **Prerequisites**: Make sure you have [napari](https://napari.org/) installed on your system. If you don't have it yet, you can easily install it using [conda-forge](https://conda-forge.org/docs/user/introduction.html):

```bash
conda create -n napari-env -c conda-forge napari pyqt python=3.12
conda activate napari-env
```

2. **Install the Package**: Next, you can install the package directly from GitHub using [pip](https://pypi.org/project/pip/):

```bash
pip install git+https://github.com/oleveque/rsar.git
```

or in developer mode:
```bash
git clone https://github.com/oleveque/rsar.git
pip install -e ./rsar
```

## Using the Python API

Here is an example of how to use the Python API to read an rSAR file and access its properties:

```python
from rsar.files import Image, Trajectory

# Load an rSAR image
image = Image('path/to/your/image.slc.toml')

# Access image properties
print(image.name)
print(image.description)
print(image.shape)
print(image.axis_units)
print(image.axis_labels)

# Load a trajectory file
trajectory = Trajectory('path/to/your/trajectory.traj.csv')

# Convert trajectory to .pos file
trajectory.to_pos_file('path/to/your/output.pos')
```

## Using the napari Plugin

1. **Launch napari**: To start using the plugin, launch napari from your terminal:

```bash
napari
```

2. **Open rSAR files**: Once napari is up and running, you can open an rSAR file by either dragging and dropping it into the napari viewer or using the command with the path to the rSAR image directly, for example:

```bash
napari <path/to/your/image.slc.toml>
```

## Issues

If you encounter any problems, please [file an issue](https://github.com/oleveque/rsar/issues) with a detailed description. Your feedback is valuable in improving the plugin and the API.