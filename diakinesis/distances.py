"""FIJI-compatible centroid to nearest six-face LMN border voxel center."""

import numpy as np
import scipy.ndimage as ndi
from scipy.spatial import cKDTree


def border6(mask):
    interior = mask.copy()
    interior[1:] &= mask[:-1]
    interior[:-1] &= mask[1:]
    interior[:,1:] &= mask[:,:-1]
    interior[:,:-1] &= mask[:,1:]
    interior[:,:,1:] &= mask[:,:,:-1]
    interior[:,:,:-1] &= mask[:,:,1:]
    interior[0] = interior[-1] = False
    interior[:,0] = interior[:,-1] = False
    interior[:,:,0] = interior[:,:,-1] = False
    return mask & ~interior


class BorderDistances:
    def __init__(self, mask, voxel_zyx_um):
        self.voxel = np.asarray(voxel_zyx_um, float)
        self.border = border6(mask)
        points = np.argwhere(self.border)
        if not len(points):
            raise ValueError('Selected LMN mask has no border voxels')
        self.tree = cKDTree(points * self.voxel)
        self.max_internal_um = float(ndi.distance_transform_edt(
            mask & ~self.border, sampling=self.voxel).max())
        self.voxel_diagonal_um = float(np.linalg.norm(self.voxel))

    def measure(self, centroid_zyx_px):
        distance = float(self.tree.query(np.asarray(centroid_zyx_px)*self.voxel)[0])
        if distance <= self.max_internal_um+self.voxel_diagonal_um:
            qc = 'pass'
        elif distance <= self.max_internal_um+2*self.voxel_diagonal_um:
            qc = 'borderline'
        else:
            qc = 'implausible'
        return distance, qc
