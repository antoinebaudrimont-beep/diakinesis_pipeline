"""Historical saved-mask and clean-stack checks; never used by production."""

import csv
import os
from pathlib import Path
import unittest

import numpy as np
import tifffile

from diakinesis.distances import BorderDistances
from diakinesis.gating import classify_objects
from diakinesis.io import inspect_dv, read_channels
from diakinesis.morphology import prepare_lmn, select_central_lmn
from diakinesis.segmentation import segment


REFERENCE=(Path(os.environ['DIAKINESIS_REFERENCE_ROOT'])
           if os.environ.get('DIAKINESIS_REFERENCE_ROOT') else None)
VALIDATION=(Path(os.environ['DIAKINESIS_VALIDATION_ROOT'])
            if os.environ.get('DIAKINESIS_VALIDATION_ROOT') else None)
FIXTURES_AVAILABLE=bool(REFERENCE and VALIDATION and
                        (REFERENCE/'to_process').is_dir() and
                        (VALIDATION/'outputs/audit_file_mapping.csv').is_file())
SAMPLES=('auxin_01','etoh_01')
VOXEL=(.2,.0643,.0643)
CHANNELS={'dapi_emission_nm':435,'lmn_emission_nm':679,'wavelength_tolerance_nm':20}
DAPI={'gaussian_radius':1,'moments_bins':20,'minimum_size_voxels':200,'fill_small_holes':False}
LMN={'gaussian_radius':1,'moments_bins':25,'minimum_size_voxels':70000}
CENTRAL={'center_gap_high_px':75,'center_gap_medium_px':35}


@unittest.skipUnless(FIXTURES_AVAILABLE,'Historical regression fixtures unavailable')
class HistoricalRegression(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with (VALIDATION/'outputs/audit_file_mapping.csv').open(newline='') as file:
            cls.mapping={r['processed_id']:r for r in csv.DictReader(file)
                         if r['status']=='matched'}

    def test_fiji_lmn_morphology_exact(self):
        for sample in SAMPLES:
            with self.subTest(sample=sample):
                labels=tifffile.imread(REFERENCE/'to_process'/f'lmn_{sample}.tif').squeeze()
                saved=tifffile.imread(REFERENCE/'to_process/binary'/f'lmn_{sample}_binary.tif').squeeze()>0
                np.testing.assert_array_equal(prepare_lmn(labels),saved)

    def test_saved_mask_cen_bor_exact(self):
        for sample in SAMPLES:
            with self.subTest(sample=sample):
                dapi=tifffile.imread(REFERENCE/'to_process'/f'dapi_{sample}.tif').squeeze()
                lmn=tifffile.imread(REFERENCE/'to_process/binary'/f'lmn_{sample}_binary.tif').squeeze()>0
                distance=BorderDistances(lmn,VOXEL)
                refs=[]
                for name in self.mapping[sample]['reference_files'].split(';'):
                    with (REFERENCE/'to_analyze'/name).open(newline='') as file:
                        refs.extend(float(row['cen-bor']) for row in csv.DictReader(file)
                                    if row['Label1'].startswith('dapi') and row['Label2']=='lmn')
                self.assertEqual(len(refs),int(dapi.max()))
                for label,expected in enumerate(refs,1):
                    center=np.argwhere(dapi==label).mean(axis=0)
                    actual,_=distance.measure(center)
                    self.assertAlmostEqual(actual,expected,delta=1e-12)

    def test_clean_stacks_retained_and_excluded(self):
        for sample,expected in (('auxin_01',6),('etoh_01',5)):
            with self.subTest(sample=sample):
                entry=self.mapping[sample]
                path=REFERENCE/entry['condition']/entry['raw_file']
                info=inspect_dv(path,CHANNELS)
                images=read_channels(info)
                lmn=segment(images['lmn'],VOXEL,LMN)
                selected,_=select_central_lmn(prepare_lmn(lmn),CENTRAL)
                dapi=segment(images['dapi'],VOXEL,DAPI)
                objects=classify_objects(dapi,selected,VOXEL,BorderDistances(selected,VOXEL))
                retained=[obj for obj in objects if obj['selected']]
                self.assertEqual(len(retained),expected)
                self.assertTrue(any(not obj['selected'] for obj in objects))


if __name__=='__main__':unittest.main()
