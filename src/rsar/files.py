from pathlib import Path
from tomlkit import load
from datetime import datetime
import numpy as np
import pandas as pd
from sargeom.coordinates import Cartographic, CartesianECEF

class Image(np.memmap):
    def __new__(cls, filename):
        filename = Path(filename)

        # Check if the specified file is a valid rSAR file
        assert filename.exists(), "File does not exist"
        assert filename.is_file(), "Path does not point to a file"
        assert filename.suffix == '.toml', "File is not a TOML file"

        # Read header information from the specified file
        with filename.open('rt', encoding='utf-8') as h:
            header = load(h)

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

        # Load the data from the file
        self = super().__new__(cls, filename.with_suffix(''),
            shape=(header['data']['row']['size'], header['data']['col']['size']),
            offset=header['data']['global_offset_byte'],
            dtype=byteorder+type,
            mode='r'
        )

        # Extract scale and offset values for data rescaling
        self._scale = self.header['data']['data']['scale']
        self._offset = self.header['data']['data']['offset']

        # Store the header information
        self._header = header

        return self
    
    def __init__(self, filename):
        pass

    def is_complex(self):
        return self._header['data']['is_complex']
    
    @property
    def axis_units(self):
        return (self._header['data']['row']['unit'], self._header['data']['col']['unit'])
    
    @property
    def axis_labels(self):
        return (self._header['data']['row']['name'], self._header['data']['col']['name'])
    
    @property
    def datetime(self):
        return datetime.fromisoformat(self._header['data']['datetime'])
    
    @property
    def name(self):
        return self._header['data']['name']
    
    @property
    def description(self):
        return self._header['data']['desc']
    
    def xlim(self):
        return (self._header['data']['row']["origin"], self._header['data']["row"]["origin"] + (self._header['data']["row"]["size"]-1) * self._header['data']["row"]["step"])
    
    def ylim(self):
        return (self._header['data']["col"]["origin"], self._header['data']["col"]["origin"] + (self._header['data']["col"]["size"]-1) * self._header['data']["col"]["step"])
    
    def xaxis(self):
        return np.linspace(*self.xlim(), self.shape[0])
    
    def yaxis(self):
        return np.linspace(*self.ylim(), self.shape[1])
    
    def is_flash(self):
        return self._header['log'].has_key('flash')
    
    def tx_integration_positions(self):
        if self.is_flash():
            start = Cartographic(
                latitude=self._header['data']['log']['flash']['azimuth']['integration']['tx']['start']['position_geo'][0],
                longitude=self._header['data']['log']['flash']['azimuth']['integration']['tx']['start']['position_geo'][1],
                altitude=self._header['data']['log']['flash']['azimuth']['integration']['tx']['start']['position_geo'][2]
            )
            center = Cartographic(
                latitude=self._header['data']['log']['flash']['azimuth']['integration']['tx']['center']['position_geo'][0],
                longitude=self._header['data']['log']['flash']['azimuth']['integration']['tx']['center']['position_geo'][1],
                altitude=self._header['data']['log']['flash']['azimuth']['integration']['tx']['center']['position_geo'][2]
            )
            end = Cartographic(
                latitude=self._header['data']['log']['flash']['azimuth']['integration']['tx']['end']['position_geo'][0],
                longitude=self._header['data']['log']['flash']['azimuth']['integration']['tx']['end']['position_geo'][1],
                altitude=self._header['data']['log']['flash']['azimuth']['integration']['tx']['end']['position_geo'][2]
            )
            return (start, center, end)
        else:
            raise ValueError("The data is not a flash integration")
    
    def tx_integration_velocities(self):
        start = CartesianECEF(
            x=self._header['data']['log']['flash']['azimuth']['integration']['tx']['start']['velocity_ecef'][0],
            y=self._header['data']['log']['flash']['azimuth']['integration']['tx']['start']['velocity_ecef'][1],
            z=self._header['data']['log']['flash']['azimuth']['integration']['tx']['start']['velocity_ecef'][2]
        )
        center = CartesianECEF(
            x=self._header['data']['log']['flash']['azimuth']['integration']['tx']['center']['velocity_ecef'][0],
            y=self._header['data']['log']['flash']['azimuth']['integration']['tx']['center']['velocity_ecef'][1],
            z=self._header['data']['log']['flash']['azimuth']['integration']['tx']['center']['velocity_ecef'][2]
        )
        end = CartesianECEF(
            x=self._header['data']['log']['flash']['azimuth']['integration']['tx']['end']['velocity_ecef'][0],
            y=self._header['data']['log']['flash']['azimuth']['integration']['tx']['end']['velocity_ecef'][1],
            z=self._header['data']['log']['flash']['azimuth']['integration']['tx']['end']['velocity_ecef'][2]
        )
        return (start, center, end)
    
    def rx_integration_positions(self):
        start = Cartographic(
            latitude=self._header['data']['log']['flash']['azimuth']['integration']['rx']['start']['position_geo'][0],
            longitude=self._header['data']['log']['flash']['azimuth']['integration']['rx']['start']['position_geo'][1],
            altitude=self._header['data']['log']['flash']['azimuth']['integration']['rx']['start']['position_geo'][2]
        )
        center = Cartographic(
            latitude=self._header['data']['log']['flash']['azimuth']['integration']['rx']['center']['position_geo'][0],
            longitude=self._header['data']['log']['flash']['azimuth']['integration']['rx']['center']['position_geo'][1],
            altitude=self._header['data']['log']['flash']['azimuth']['integration']['rx']['center']['position_geo'][2]
        )
        end = Cartographic(
            latitude=self._header['data']['log']['flash']['azimuth']['integration']['rx']['end']['position_geo'][0],
            longitude=self._header['data']['log']['flash']['azimuth']['integration']['rx']['end']['position_geo'][1],
            altitude=self._header['data']['log']['flash']['azimuth']['integration']['rx']['end']['position_geo'][2]
        )
        return (start, center, end)
    
    def rx_integration_velocities(self):
        start = CartesianECEF(
            x=self._header['data']['log']['flash']['azimuth']['integration']['rx']['start']['velocity_ecef'][0],
            y=self._header['data']['log']['flash']['azimuth']['integration']['rx']['start']['velocity_ecef'][1],
            z=self._header['data']['log']['flash']['azimuth']['integration']['rx']['start']['velocity_ecef'][2]
        )
        center = CartesianECEF(
            x=self._header['data']['log']['flash']['azimuth']['integration']['rx']['center']['velocity_ecef'][0],
            y=self._header['data']['log']['flash']['azimuth']['integration']['rx']['center']['velocity_ecef'][1],
            z=self._header['data']['log']['flash']['azimuth']['integration']['rx']['center']['velocity_ecef'][2]
        )
        end = CartesianECEF(
            x=self._header['data']['log']['flash']['azimuth']['integration']['rx']['end']['velocity_ecef'][0],
            y=self._header['data']['log']['flash']['azimuth']['integration']['rx']['end']['velocity_ecef'][1],
            z=self._header['data']['log']['flash']['azimuth']['integration']['rx']['end']['velocity_ecef'][2]
        )
        return (start, center, end)

class Trajectory(pd.DataFrame):
    def __init__(self, filename):
        filename = Path(filename)

        # Check if the specified file is a valid trajectory file
        assert filename.exists(), "File does not exist"
        assert filename.is_file(), "Path does not point to a file"
        assert filename.suffix == '.traj.csv', "File is not a .traj.csv file"

        # Read trajectory data from the specified file
        data = pd.read_csv(filename, delimiter=';', comment='#')

        # Check if the file has the expected columns
        expected_columns = ['TIMESTAMP_S', 'LON_WGS84_DEG', 'LAT_WGS84_DEG', 'HEIGHT_WGS84_M', 'HEADING_DEG', 'ELEVATION_DEG', 'BANK_DEG']
        if not all(column in data.columns for column in expected_columns):
            raise ValueError("The .traj.csv file must contain the columns TIMESTAMP_S, LON_WGS84_DEG, LAT_WGS84_DEG, HEIGHT_WGS84_M, HEADING_DEG, ELEVATION_DEG, BANK_DEG")

        super().__init__(data)