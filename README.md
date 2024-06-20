# napari-rsar

A plugin for napari that facilitates the visualization and analysis of rSAR files.

## Installation

1. **Prerequisites**: Make sure you have [napari](https://napari.org/) installed on your system. If not, you can install it using [conda-forge](https://conda-forge.org/docs/user/introduction.html):

```bash
conda create -n napari-env -c conda-forge napari pyqt python=3.11
conda activate napari-env
```

2. **Install the Plugin**: You can install the plugin directly from GitHub using [pip](https://pypi.org/project/pip/):

```bash
pip install git+https://github.com/oleveque/napari-rsar.git
```

or in developer mode:
```bash
git clone https://github.com/oleveque/napari-rsar.git
pip install -e ./napari-rsar
```

3. **Launch napari**: Launch napari from your terminal:

```bash
napari
```

4. **Open rSAR files**: Once napari is launched, you can open a rSAR file by dragging and dropping it to the napari viewer.

## Issues

If you encounter any problems, please [file an issue](https://github.com/oleveque/napari-rsar/issues) with a detailed description. Your feedback is valuable in improving the plugin.