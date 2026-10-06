"""Retain central-nucleus DAPI and trace every full-field detection."""

import math

import numpy as np
import scipy.ndimage as ndi


def classify_objects(labels, selected_lmn, voxel_zyx_um, distances):
    counts = np.bincount(labels.ravel())
    ids = np.flatnonzero(counts)
    ids = ids[ids > 0]
    if not len(ids):
        raise ValueError('DAPI segmentation found no objects')
    coords = np.nonzero(labels)
    values = labels[coords]
    sums = [np.bincount(values, weights=axis, minlength=len(counts)) for axis in coords]
    inside_counts = np.bincount(labels.ravel(), weights=selected_lmn.ravel(),
                                 minlength=len(counts))
    outside = ndi.distance_transform_edt(~selected_lmn, sampling=voxel_zyx_um)
    nearest = ndi.minimum(outside, labels, index=ids)
    shape = np.asarray(labels.shape)
    voxel_volume = math.prod(voxel_zyx_um)
    diagonal = math.sqrt(sum(v*v for v in voxel_zyx_um))
    result=[]
    for label, minimum in zip(ids, nearest):
        volume = int(counts[label])
        center = np.asarray([s[label]/volume for s in sums])
        center_index = tuple(np.clip(np.rint(center).astype(int),0,shape-1))
        center_inside = bool(selected_lmn[center_index])
        center_outside = float(outside[center_index])
        volume_um3 = volume*voxel_volume
        radius = (3*volume_um3/(4*math.pi))**(1/3)
        borderline = (not center_inside and float(minimum)<=diagonal and
                      center_outside<=radius+diagonal)
        classification = ('inside' if center_inside else
                          'borderline_retained' if borderline else 'outside_excluded')
        distance,geometry_qc=distances.measure(center)
        reason = ('centroid_inside_selected_lmn' if center_inside else
                  'touches_selected_lmn_with_centroid_near_border' if borderline else
                  'outside_selected_lmn_tolerance')
        result.append(dict(dapi_object_id=int(label), selected=classification!='outside_excluded',
            selection_class=classification, selection_reason=reason,
            volume_voxels=volume, volume_um3=volume_um3,
            centroid_z_px=float(center[0]), centroid_y_px=float(center[1]),
            centroid_x_px=float(center[2]), centroid_inside_lmn=center_inside,
            fraction_inside_lmn=float(inside_counts[label]/volume),
            minimum_distance_to_lmn_um=float(minimum),
            centroid_outside_lmn_um=center_outside,
            cen_bor_um=distance, geometry_qc=geometry_qc))
    return result
