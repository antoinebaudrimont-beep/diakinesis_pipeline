"""Stable CSV schemas and mask writing."""

import csv

import tifffile


FIELDS = ['source_file','condition','nucleus_id','selected_lmn_component',
          'dapi_object_id','auto_selected','selected','manual_selected','manual_comment',
          'selection_class','selection_reason','histogram_bins','volume_voxels','volume_um3',
          'centroid_x','centroid_y','centroid_z',
          'centroid_x_px','centroid_y_px','centroid_z_px','centroid_inside_lmn',
          'fraction_inside_lmn','minimum_distance_to_lmn_um','centroid_outside_lmn_um',
          'cen_bor_um','geometry_qc','review_status','voxel_size_x_um',
          'voxel_size_y_um','voxel_size_z_um']


def write_csv(path,rows,fields):
    with path.open('w',newline='') as file:
        writer=csv.DictWriter(file,fieldnames=fields,extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def save_masks(directory,sample,labels,lmn_labels,prepared):
    tifffile.imwrite(directory/f'{sample}_dapi_labels.tif',labels.astype('int32'),
                     photometric='minisblack')
    tifffile.imwrite(directory/f'{sample}_lmn_mask.tif',(lmn_labels>0).astype('uint8'),
                     photometric='minisblack')
    tifffile.imwrite(directory/f'{sample}_lmn_prepared.tif',prepared.astype('uint8'),
                     photometric='minisblack')
