"""Whole-stack PartSeg/SimpleITK-equivalent Moments segmentation."""

import numpy as np
import SimpleITK as sitk


def segment(image, voxel_zyx_um, settings):
    if settings.get('fill_small_holes', False):
        raise ValueError('Small-hole filling is not implemented in the validated profile')
    radius = float(settings['gaussian_radius'])
    bins = int(settings['moments_bins'])
    minimum = int(settings['minimum_size_voxels'])
    if radius < 0 or bins < 2 or minimum < 1:
        raise ValueError('Invalid segmentation parameters')
    z,y,x = voxel_zyx_um
    base = min(z,y,x)
    variance_xyz = [radius*base/x, radius*base/y, radius*base/z]
    filtered = sitk.DiscreteGaussian(sitk.GetImageFromArray(image), variance_xyz)
    binary = sitk.MomentsThreshold(filtered, 0, 1, bins)
    connected = sitk.ConnectedComponent(binary, True)
    labels = sitk.GetArrayFromImage(sitk.RelabelComponent(connected, 20))
    counts = np.bincount(labels.ravel())
    keep = np.flatnonzero(counts >= minimum)
    keep = keep[keep > 0]
    return np.where(np.isin(labels, keep), labels, 0).astype(np.int32)
