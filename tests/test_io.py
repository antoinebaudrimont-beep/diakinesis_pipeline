"""Raw DV metadata must resolve channels and voxel sizes without reference files."""

from pathlib import Path
import struct
import tempfile
import unittest

from diakinesis.io import discover, inspect_dv
from diakinesis.qc import distinct_colors


CHANNELS={'dapi_emission_nm':435,'lmn_emission_nm':679,'wavelength_tolerance_nm':20}


class RawInputChecks(unittest.TestCase):
    def make_dv(self,folder,wavelengths=(679,435)):
        path=Path(folder)/'synthetic_R3D_D3D.dv'
        header=bytearray(1024)
        struct.pack_into('<4i',header,0,2,2,4,2)
        struct.pack_into('<3f',header,40,.064346,.064346,.2)
        struct.pack_into('<i',header,92,0)
        struct.pack_into('<h',header,180,1)
        struct.pack_into('<h',header,182,2)
        struct.pack_into('<h',header,196,2)
        struct.pack_into('<2h',header,198,*wavelengths)
        path.write_bytes(header+bytearray(4*2*2*4))
        return path

    def test_reads_embedded_channels_and_calibration(self):
        with tempfile.TemporaryDirectory() as folder:
            path=self.make_dv(folder)
            self.assertEqual(discover(folder),[path])
            info=inspect_dv(path,CHANNELS)
            self.assertEqual(info['channel_indices'],{'dapi':1,'lmn':0})
            self.assertAlmostEqual(info['metadata_voxel_size_zyx_um'][0],.2,places=6)

    def test_rejects_missing_channel_metadata(self):
        with tempfile.TemporaryDirectory() as folder:
            path=self.make_dv(folder,wavelengths=(679,0))
            with self.assertRaisesRegex(ValueError,'Missing channel wavelengths'):
                inspect_dv(path,CHANNELS)

    def test_rejects_empty_folder(self):
        with tempfile.TemporaryDirectory() as folder:
            with self.assertRaisesRegex(ValueError,'No .* files found'):
                discover(folder)

    def test_selected_qc_colors_remain_distinct_for_large_counts(self):
        colors=distinct_colors(30)
        self.assertEqual(len(colors),len(set(colors)))
        self.assertNotIn('#e31a1c',colors)


if __name__=='__main__':unittest.main()
