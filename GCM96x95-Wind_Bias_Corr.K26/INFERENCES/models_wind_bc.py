import torch
import numpy as np
from pathlib import Path
from unet3d_inference import UNet3D


# ++++++++++++++++++++++++++++ #
#   USER-DEFINED PARAMETERS    #
# ++++++++++++++++++++++++++++ #
# Normalisation and checkpoint files
NORM_FILE_GLOB = '*_normalisation.npz'
CHECKPOINT_FILE_GLOB = '*_state_dict.pt'

# Normalisation mode for inputs:
#   'pointwise' : z-score grid point
#   'level'     : mean/std profiles by level
#   'global'    : unique scalars
#   'none'      : no normalisation -- ONLY option to use a different grid
NORM_MODE = 'pointwise'

# Topography normalisation scalars
# If None, field is self-normalized ('global'-type)
TOPO_MEAN = None
TOPO_STD = None

# Torch dtype used for both inputs and outputs.
DATA_DTYPE = 'float32'

# ++++++++++++++++++++++++++++ #
#         Module state          #
# ++++++++++++++++++++++++++++ #
# Channel order used at training time
FEATURE_ORDER = ['pres', 'psol', 'tpot', 'u', 'v', 'topo', 'cos_day', 'sin_day']
POINTWISE_3D_CHANNELS = ('pres', 'u', 'v', 'tpot')   
POINTWISE_2D_CHANNELS = ('psol',)                    

if NORM_MODE not in ('pointwise', 'level', 'global', 'none'):
    raise ValueError(f"NORM_MODE invalide: {NORM_MODE!r} (pointwise|level|global|none)")

_model: torch.nn.Module | None = None
_norm: dict | None = None


def _load_normalisation() -> dict | None:
    files = sorted(Path('.').glob(NORM_FILE_GLOB))
    if not files:
        print(f"[unet_inference_online] WARNING: {NORM_FILE_GLOB} not found ")
        return None
    data = np.load(files[0])
    norm = {key: data[key] for key in data.files}
    print(f'[unet_inference_online] Normalisation loaded from: {files[0]}')
    return norm


def _load_state_dict_checkpoint(checkpoint_path: Path) -> tuple[torch.nn.Module, dict]:
    try:
        payload = torch.load(checkpoint_path, map_location='cpu', weights_only=True)
    except Exception:
        payload = torch.load(checkpoint_path, map_location='cpu', weights_only=False)
    config = payload.get('config')
    if config is None or 'state_dict' not in payload:
        raise ValueError(f'not a converted state dict file: {checkpoint_path}')
    model = UNet3D(attention_gate=config['attention_gate'],
                   in_channels=config['in_channels'],
                   out_channels=config['out_channels'],
                   base_channels=config['base_channels'])
    model.load_state_dict(payload['state_dict'], strict=True)
    return model, payload.get('meta', {})


def init_model(checkpoint_path: Path | None = None) -> torch.nn.Module:
    global _model
    if _model is not None:
        return _model
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    if checkpoint_path is None:
        ckpts = sorted(Path('.').glob(CHECKPOINT_FILE_GLOB))
        checkpoint_path = ckpts[0] if ckpts else None
    if checkpoint_path is not None and checkpoint_path.exists():
        print(f'[unet_inference_online] Loading trained model from: {checkpoint_path}')
        model, _meta = _load_state_dict_checkpoint(checkpoint_path)
    else:
        print(f"[unet_inference_online] WARNING: {CHECKPOINT_FILE_GLOB} not found ")
        model = UNet3D(attention_gate=False, in_channels=len(FEATURE_ORDER), out_channels=2)
    for param in model.parameters():
        param.requires_grad = False
    model.eval()
    model.to(device)
    _model = model
    return _model


init_model()
_norm = _load_normalisation()


# ------------------------------ #
#       Preprocessing utils      #
# ------------------------------ #
def _zscore(field: np.ndarray, mean, std) -> np.ndarray:
    std = np.asarray(std)
    std_safe = np.where(std > 0, std, 1.0)
    out = np.where(std > 0, (field - mean) / std_safe, 0.0)
    return out.astype(DATA_DTYPE)


def _pad_axis_duplicate_last(array: np.ndarray, axis: int, multiple: int = 8) -> np.ndarray:
    size = array.shape[axis]
    pad = (-size) % multiple
    if pad == 0:
        return array
    slicer = [slice(None)] * array.ndim
    slicer[axis] = slice(-1, None)
    return np.concatenate([array, np.repeat(array[tuple(slicer)], pad, axis=axis)], axis=axis)


def _to_lvl_lat_lon(field) -> np.ndarray:
    return np.asarray(field).transpose(2, 1, 0)[None].astype(DATA_DTYPE)


def _stats(varname: str, mode: str):
    m = _norm[f'{varname}_mean'].astype(np.float64)
    s = _norm[f'{varname}_std'].astype(np.float64)
    if mode == 'level' and varname != 'psol':
        var = (s ** 2).mean(axis=(1, 2)) + m.var(axis=(1, 2))
        return m.mean(axis=(1, 2)).astype(np.float32), np.sqrt(var).astype(np.float32)
    var = (s ** 2).mean() + m.var()
    return np.float32(m.mean()), np.float32(np.sqrt(var))


def _apply_normalisation(fields: dict) -> None:
    lvl, lat, lon = fields['u'].shape[1:]
    nsrc, nlat, nlon = _norm['u_mean'].shape
    if (lat, lon) != (nlat, nlon) or lvl not in (nsrc, nsrc - 1):
        raise ValueError(
            f'Grid ({lvl},{lat},{lon}) incompatible with normalisation '
            f'(Expected : {nsrc - 1} or {nsrc} levels x {nlat} x {nlon}). '
            "Use NORM_MODE = 'none' for another grid (physical consistency not guaranteed).")
    sl = slice(None) if lvl == nsrc else slice(None, -1)
    if NORM_MODE == 'pointwise':
        for v in POINTWISE_3D_CHANNELS:
            fields[v] = _zscore(fields[v], _norm[f'{v}_mean'][sl], _norm[f'{v}_std'][sl])
        fields['psol'] = _zscore(fields['psol'], _norm['psol_mean'], _norm['psol_std'])
    elif NORM_MODE == 'level':
        for v in POINTWISE_3D_CHANNELS:
            m, s = _stats(v, 'level')                     
            fields[v] = _zscore(fields[v],
                                m[:lvl][None, :, None, None], s[:lvl][None, :, None, None])
        m, s = _stats('psol', 'level')                    
        fields['psol'] = _zscore(fields['psol'], m, s)
    else:  # 'global'
        for v in POINTWISE_3D_CHANNELS + POINTWISE_2D_CHANNELS:
            m, s = _stats(v, 'global')                    
            fields[v] = _zscore(fields[v], m, s)
    # topo: user scalars if provided, else 'global'-type self-normalization
    if TOPO_MEAN is not None and TOPO_STD is not None:
        topo_mean, topo_std = TOPO_MEAN, TOPO_STD
    else:
        topo_mean = fields['topo'][0, 0].mean()
        topo_std = fields['topo'][0, 0].std()
    fields['topo'] = _zscore(fields['topo'], topo_mean, topo_std)


# ------------------------------ #
#       Input Formatting         #
# ------------------------------ #
def Is_None(*inputs):
    """ Test presence of at least one None in inputs """
    return any(item is None for item in inputs)


def format_inputs(pres, u, v, tpot, psol, topo, cos_day, sin_day):
    """
    Reproduit le preprocessing d'entraînement puis empile les canaux.

    Entrees (convention physique LMDZ, un pas de temps) : (Lon, Lat, Levels) ;
    psol et topo ont une dimension de niveau singleton : (Lon, Lat, 1).

    Pipeline : permutation -> normalisation selon NORM_MODE ('none' : aucune)
    -> padding des niveaux vers un multiple de 8 (79 -> 80) -> empilement
    channels-last (1, Levels, Lat, Lon, 8).
    """
    fields = {
        'pres': _to_lvl_lat_lon(pres),
        'u': _to_lvl_lat_lon(u),
        'v': _to_lvl_lat_lon(v),
        'tpot': _to_lvl_lat_lon(tpot),
        'psol': _to_lvl_lat_lon(psol),      # (1, 1, Lat, Lon)
        'topo': _to_lvl_lat_lon(topo),      # (1, 1, Lat, Lon)
        'cos_day': _to_lvl_lat_lon(cos_day),
        'sin_day': _to_lvl_lat_lon(sin_day),
    }

    if _norm is not None and NORM_MODE != 'none':
        _apply_normalisation(fields)

    # Padding des niveaux vers un multiple de 8 (duplication du pseudo-niveau).
    for name, field in fields.items():
        if name not in ('psol', 'topo'):
            fields[name] = _pad_axis_duplicate_last(field, axis=1)

    # Broadcast vers la forme commune puis empilement channels-last.
    ref_shape = fields['u'].shape
    arrs = [np.broadcast_to(fields[name], ref_shape) for name in FEATURE_ORDER]
    return np.stack(arrs, axis=-1).astype(DATA_DTYPE)


# ------------------------------ #
#       Main Model Routine       #
# ------------------------------ #
def inference_pytorch_lightning(u, v, tpot, pres, psol, topo, cos_day, sin_day):
    """
    Compute air-sea momentum corrections (du, dv) from the trained ANN.

    Single time-step online inference. All numpy inputs use the
    (Lon, Lat, Levels) axis convention (LMDZ physics).
    psol and topo have a singleton level dimension: (Lon, Lat, 1).

    Les entrees sont normalisees selon NORM_MODE et paddees en niveau ; les
    sorties sont recadrees sur la grille d'entree. Aucune denormalisation en
    sortie : du/dv sont deja en unites physiques (increment de nudging en vent
    naturel).

    Returns:
      du, dv : each (Lon, Lat, Levels) as float32 numpy arrays.

    """
    if Is_None(u, v, tpot, pres, psol, topo, cos_day, sin_day):
        return None, None

    model = _model
    lon_size, lat_size, lvl_size = int(u.shape[0]), int(u.shape[1]), int(u.shape[2])

    stacked = format_inputs(pres, u, v, tpot, psol, topo, cos_day, sin_day)
    inputs = torch.from_numpy(stacked.copy()).permute(0, 4, 1, 2, 3).contiguous()

    # Predict -- (1, 2, Levels_p, Lat_p, Lon_p)
    with torch.no_grad():
        preds = model(inputs).permute(0, 2, 3, 4, 1).cpu().numpy()

    # (Levels_p, Lat_p, Lon_p) -> (Lon_p, Lat_p, Levels_p), recadrage du padding.
    du = preds[0, ..., 0].transpose(2, 1, 0)[:lon_size, :lat_size, :lvl_size].astype(DATA_DTYPE)
    dv = preds[0, ..., 1].transpose(2, 1, 0)[:lon_size, :lat_size, :lvl_size].astype(DATA_DTYPE)
    return np.ascontiguousarray(du), np.ascontiguousarray(dv)



if __name__ == '__main__':

    # Test with dummy inputs 
    nsrc, lat, lon = _norm['u_mean'].shape
    raw = {v: np.ascontiguousarray(_norm[f'{v}_mean'].transpose(2, 1, 0)) for v in POINTWISE_3D_CHANNELS}                               
    raw['psol'] = np.ascontiguousarray(_norm['psol_mean'].T[..., None])      
    raw['topo'] = np.zeros((lon, lat, 1), dtype=DATA_DTYPE)                  
    raw['cos_day'] = np.full((lon, lat, nsrc), 0.7 * np.sqrt(2.0), dtype=DATA_DTYPE)
    raw['sin_day'] = np.full((lon, lat, nsrc), 0.5 * np.sqrt(2.0), dtype=DATA_DTYPE)
    j0, j1 = lat // 4, lat - lat // 4
    raw['u'][:, j0:j1, 0] += 3.0                                              
    raw['tpot'][:, j0:j1, 0] += 2.0

    # Call NN
    du, dv = inference_pytorch_lightning(raw['u'], raw['v'], raw['tpot'], raw['pres'],
                                         raw['psol'], raw['topo'], raw['cos_day'], raw['sin_day'])

    # Check results
    assert du.shape == (lon, lat, nsrc), f'du shape mismatch: {du.shape}'
    print(f'du [{du.min():+.2e},{du.max():+.2e}] mean {du.mean():+.2e} '
          f'|surf|max {np.abs(du[..., 0]).max():.2e}')
    print(f'dv [{dv.min():+.2e},{dv.max():+.2e}] mean {dv.mean():+.2e} '
          f'|surf|max {np.abs(dv[..., 0]).max():.2e}')
    print('Test successful')
