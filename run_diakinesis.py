#!/usr/bin/env python3
"""Analyze a folder of raw two-channel deconvolved DeltaVision DV stacks."""

import argparse
import csv
from datetime import datetime, timezone
import importlib.metadata
import logging
from pathlib import Path
import sys

import numpy as np
import yaml

from diakinesis import __version__
from diakinesis.io import discover, inspect_dv, read_channels
from diakinesis.segmentation import segment
from diakinesis.morphology import prepare_lmn, select_central_lmn
from diakinesis.distances import BorderDistances
from diakinesis.gating import classify_objects
from diakinesis.qc import save_qc
from diakinesis.output import FIELDS, save_masks, write_csv


def parse_args(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,help='Folder containing *_R3D_D3D.dv')
    parser.add_argument('--output',type=Path,help='New results folder (default: input/diakinesis_results)')
    parser.add_argument('--config',type=Path,default=Path(__file__).with_name('config.yaml'))
    parser.add_argument('--condition',help='Condition label (default: input folder name)')
    parser.add_argument('--dapi-bins',type=int,help='Whole-stack Moments bins for every input file')
    parser.add_argument('--calibration',choices=['metadata','legacy'],help='Voxel calibration mode')
    return parser.parse_args(argv)


def choose_input_folder():
    try:
        import tkinter as tk
        from tkinter import filedialog
        root=tk.Tk();root.withdraw()
        name=filedialog.askdirectory(title='Select folder containing deconvolved DV stacks')
        root.destroy()
    except Exception as exc:
        raise ValueError('Folder chooser could not open; supply --input on the command line') from exc
    if not name:
        raise ValueError('No input folder selected')
    return Path(name)


def versions():
    names=['numpy','scipy','SimpleITK','tifffile','matplotlib','PyYAML']
    return {name:importlib.metadata.version(name) for name in names}


def run(argv=None):
    args=parse_args(argv)
    source=(args.input or choose_input_folder()).expanduser().resolve()
    if not source.is_dir():
        raise ValueError(f'Input folder does not exist: {source}')
    config=yaml.safe_load(args.config.read_text())
    if args.dapi_bins is not None:
        config['dapi']['moments_bins']=args.dapi_bins
    if args.calibration:
        config['calibration']['mode']=args.calibration
    if config['calibration']['mode'] not in ('metadata','legacy'):
        raise ValueError('calibration.mode must be metadata or legacy')
    if config['gating']['mode']!='tolerant':
        raise ValueError('Only validated tolerant gating is supported')
    condition=args.condition or config.get('condition') or source.name
    config['condition']=condition
    files=discover(source)
    info=[inspect_dv(path,config['channels']) for path in files]
    dest=(args.output or source/'diakinesis_results').expanduser().resolve()
    if dest == source or dest in (p.resolve() for p in files):
        raise ValueError('Output must be a separate results folder')
    if dest.exists() and any(dest.iterdir()):
        raise ValueError(f'Output folder already contains files: {dest}. Choose a new --output.')
    dest.mkdir(parents=True,exist_ok=True)
    for sub in ('qc','masks','tables','logs'):
        (dest/sub).mkdir(exist_ok=True)
    log_path=dest/'analysis_log.txt'
    logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s',
                        handlers=[logging.FileHandler(log_path),logging.StreamHandler(sys.stdout)],force=True)
    log=logging.getLogger('diakinesis')
    log.info('Pipeline version %s; started %s',__version__,datetime.now(timezone.utc).isoformat())
    log.info('Input %s; output %s; condition %s',source,dest,condition)
    log.info('Found %d DV files: %s',len(files),', '.join(p.name for p in files))
    log.info('Software versions: %s',versions())
    (dest/'logs/environment.txt').write_text(
        f'Pipeline version: {__version__}\nPython: {sys.version}\n'+
        ''.join(f'{name}: {version}\n' for name,version in versions().items()))
    (dest/'config_used.yaml').write_text(yaml.safe_dump(config,sort_keys=False))
    inventory=[]
    for item in info:
        raw=Path(item['path']).name
        metadata=item['metadata_voxel_size_zyx_um']
        used=(metadata if config['calibration']['mode']=='metadata' else
              tuple(float(config['calibration']['legacy_voxel_size_um'][axis]) for axis in 'zyx'))
        item['used_voxel_size_zyx_um']=used
        log.info('%s channels=%s nm roles=%s metadata voxel z/y/x=%s µm used=%s µm',
                 raw,item['wavelengths_nm'],item['channel_indices'],metadata,used)
        inventory.append(dict(source_file=raw,channel_wavelengths_nm=';'.join(map(str,item['wavelengths_nm'])),
                              dapi_channel_index=item['channel_indices']['dapi'],
                              lmn_channel_index=item['channel_indices']['lmn'],
                              metadata_z_um=metadata[0],metadata_y_um=metadata[1],metadata_x_um=metadata[2],
                              used_z_um=used[0],used_y_um=used[1],used_x_um=used[2],
                              calibration_mode=config['calibration']['mode']))
    write_csv(dest/'tables/input_inventory.csv',inventory,list(inventory[0]))
    all_rows=[];stacks=[]
    for number,item in enumerate(info,1):
        raw=Path(item['path']).name
        log.info('Processing %d/%d %s',number,len(info),raw)
        try:
            channels=read_channels(item)
            voxel=item['used_voxel_size_zyx_um']
            lmn_labels=segment(channels['lmn'],voxel,config['lmn'])
            prepared=prepare_lmn(lmn_labels,**dict(
                dilation_iterations=config['morphology']['dilation_iterations'],
                hole_fill_iterations=config['morphology']['hole_fill_iterations']))
            selected_lmn,lmn=select_central_lmn(prepared,config['lmn'])
            dapi_labels=segment(channels['dapi'],voxel,config['dapi'])
            distances=BorderDistances(selected_lmn,voxel)
            objects=classify_objects(dapi_labels,selected_lmn,voxel,distances)
            retained=[obj for obj in objects if obj['selected']]
            excluded=[obj for obj in objects if not obj['selected']]
            review=[]
            if lmn['selection_confidence']!='high':review.append('ambiguous_lmn_component')
            if len(retained)<int(config['review']['low_count']):review.append('few_retained_dapi')
            if len(retained)>int(config['review']['high_count']):review.append('many_retained_dapi')
            if any(obj['selection_class']=='borderline_retained' for obj in retained):
                review.append('borderline_dapi_gate')
            if any(obj['geometry_qc']!='pass' for obj in retained):
                review.append('geometry_qc')
            if len(retained)>2:
                volumes=np.asarray([obj['volume_voxels'] for obj in retained])
                if max(volumes)>float(config['review']['merge_volume_ratio'])*np.median(volumes):
                    review.append('possible_merged_dapi')
            status='needs_review' if review else 'auto_pass'
            for obj in objects:
                all_rows.append(dict(source_file=raw,condition=condition,
                    nucleus_id=lmn['selected_lmn_component'],
                    selected_lmn_component=lmn['selected_lmn_component'],
                    auto_selected=obj['selected'],manual_selected='',manual_comment='',
                    histogram_bins=config['dapi']['moments_bins'],review_status=status,
                    voxel_size_z_um=voxel[0],voxel_size_y_um=voxel[1],voxel_size_x_um=voxel[2],**obj))
                all_rows[-1].update(centroid_x=obj['centroid_x_px'],
                                    centroid_y=obj['centroid_y_px'],centroid_z=obj['centroid_z_px'])
            sample=raw.removesuffix('_R3D_D3D.dv')
            if config['qc']['save_masks']:
                save_masks(dest/'masks',sample,dapi_labels,lmn_labels,selected_lmn)
            if config['qc']['save_png']:
                save_qc(dest/'qc'/f'{sample}_qc.png',raw,channels['dapi'],dapi_labels,
                        selected_lmn,objects,config['dapi']['moments_bins'],lmn['selection_confidence'])
            stacks.append(dict(source_file=raw,condition=condition,selected_lmn_component=lmn['selected_lmn_component'],
                lmn_components_detected=lmn['lmn_components_detected'],
                lmn_selection_confidence=lmn['selection_confidence'],retained=len(retained),excluded=len(excluded),
                review_status=status,review_reasons=';'.join(review),
                maximum_internal_distance_um=distances.max_internal_um))
            log.info('%s selected LMN=%s confidence=%s retained=%d excluded=%d review=%s',
                     raw,lmn['selected_lmn_component'],lmn['selection_confidence'],len(retained),len(excluded),
                     ';'.join(review) or 'none')
            if review:log.warning('%s requires visual review: %s',raw,';'.join(review))
        except Exception:
            log.exception('Processing failed for %s',raw)
            raise
    write_csv(dest/'measurements_all.csv',all_rows,FIELDS)
    chosen=[row for row in all_rows if row['selected']]
    write_csv(dest/'measurements_selected.csv',chosen,FIELDS)
    compatibility=[{'filename':row['source_file'],'Label1':f"dapi{row['dapi_object_id']}",
                    'Label2':'lmn','cen-bor':row['cen_bor_um']} for row in chosen]
    write_csv(dest/'summary_long.csv',compatibility,['filename','Label1','Label2','cen-bor'])
    write_csv(dest/'tables/stack_qc.csv',stacks,list(stacks[0]))
    log.info('Finished %d stacks; %d detected DAPI, %d retained, %d excluded',
             len(stacks),len(all_rows),len(chosen),len(all_rows)-len(chosen))
    return dest


if __name__=='__main__':
    try:run()
    except (ValueError,KeyError,OSError) as exc:
        print(f'Error: {exc}',file=sys.stderr)
        sys.exit(2)
