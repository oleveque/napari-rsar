# napari-rsar

A plugin for napari that facilitates the visualization and analysis of rSAR files.

## Installation

1. **Prerequisites**: Make sure you have [napari](https://napari.org/) installed on your system. If you don't have it yet, you can easily install it using [conda-forge](https://conda-forge.org/docs/user/introduction.html):

```bash
conda create -n napari-env -c conda-forge napari pyqt python=3.11
conda activate napari-env
```

2. **Install the Plugin**: Next, you can install the plugin directly from GitHub using [pip](https://pypi.org/project/pip/):

```bash
pip install git+https://github.com/oleveque/napari-rsar.git
```

or in developer mode:
```bash
git clone https://github.com/oleveque/napari-rsar.git
pip install -e ./napari-rsar
```

3. **Launch napari**: To start using the plugin, launch napari from your terminal:

```bash
napari
```

4. **Open rSAR files**: Once napari is up and running, you can open an rSAR file by either dragging and dropping it into the napari viewer or using the command with the path to the rSAR image directly, for example:

```bash
napari <path/to/your/image.slc.toml>
```

## Issues

If you encounter any problems, please [file an issue](https://github.com/oleveque/napari-rsar/issues) with a detailed description. Your feedback is valuable in improving the plugin.