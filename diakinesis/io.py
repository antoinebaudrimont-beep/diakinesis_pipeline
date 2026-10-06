"""Read the audited two-channel DeltaVision float32 DV layout and metadata."""

from pathlib import Path
import struct

import numpy as np


def discover(folder):
    files = sorted(Path(folder).glob('*_R3D_D3D.dv'))
    if not files:
        raise ValueError(f'No *_R3D_D3D.dv files found in {folder}')
    return files


def inspect_dv(path, channel_config):
    path = Path(path)
    with path.open('rb') as file:
        header = file.read(1024)
    if len(header) != 1024:
        raise ValueError(f'Incomplete DV header: {path.name}')
    nx, ny, sections, mode = struct.unpack_from('<4i', header)
    sx, sy, sz = struct.unpack_from('<3f', header, 40)
    extended = struct.unpack_from('<i', header, 92)[0]
    times, sequence, waves = (struct.unpack_from('<h', header, offset)[0]
                              for offset in (180, 182, 196))
    if mode != 2 or min(nx, ny, sections) <= 0 or extended < 0:
        raise ValueError(f'Unsupported DV layout or pixel type: {path.name}')
    if times != 1 or sequence != 2 or waves != 2 or sections % waves:
        raise ValueError(f'Expected one time point and two channel Z blocks: {path.name}; '
                         f'found times={times}, sequence={sequence}, channels={waves}')
    if not all(np.isfinite([sx, sy, sz])) or min(sx, sy, sz) <= 0:
        raise ValueError(f'Missing or invalid physical voxel size: {path.name}')
    wavelengths = [struct.unpack_from('<h', header, 198 + 2*i)[0]
                   for i in range(waves)]
    if any(w <= 0 for w in wavelengths):
        raise ValueError(f'Missing channel wavelengths in DV metadata: {path.name}')
    tolerance = float(channel_config['wavelength_tolerance_nm'])
    roles = {}
    for role in ('dapi', 'lmn'):
        expected = float(channel_config[f'{role}_emission_nm'])
        matches = [i for i, w in enumerate(wavelengths) if abs(w-expected) <= tolerance]
        if len(matches) != 1:
            raise ValueError(f'{path.name}: {role} channel ambiguous/missing; '
                             f'DV wavelengths={wavelengths} nm, expected {expected}±{tolerance} nm. '
                             'Set channels in config.yaml; do not guess channel order.')
        roles[role] = matches[0]
    if roles['dapi'] == roles['lmn']:
        raise ValueError(f'{path.name}: DAPI and LMN resolve to the same channel')
    expected_bytes = 1024 + extended + sections*ny*nx*4
    if path.stat().st_size < expected_bytes:
        raise ValueError(f'Incomplete DV image data: {path.name}')
    return dict(path=str(path), shape_zyx=(sections//waves, ny, nx),
                extension_bytes=extended, wavelengths_nm=wavelengths,
                channel_indices=roles,
                metadata_voxel_size_zyx_um=(float(sz), float(sy), float(sx)))


def read_channels(info):
    z, ny, nx = info['shape_zyx']
    data = np.memmap(info['path'], dtype='<f4', mode='r',
                     offset=1024+info['extension_bytes'], shape=(2, z, ny, nx))
    return {role:data[index] for role,index in info['channel_indices'].items()}
