import os

from .lib.sutils import *

ROOT_DIR = Path(__file__).resolve().parent
LIB_DIR = ROOT_DIR / 'lib'
DATA_DIR = ROOT_DIR.parent / 'data'

# if mosek env var is not given check custom valid paths
_MOSEK_LIC_FILE_ENV = 'MOSEKLM_LICENSE_FILE'
MOSEK_SUPPORTED_PATHS = ['~/mosek/mosek.lic', '~/.mosek/mosek.lic', '~/.mosek.lic']
if _MOSEK_LIC_FILE_ENV not in os.environ:
    for p in MOSEK_SUPPORTED_PATHS:
        _p = Path(p).expanduser()
        if _p.exists():
            _p = str(_p)
            os.environ[_MOSEK_LIC_FILE_ENV] = _p
            os.putenv(_MOSEK_LIC_FILE_ENV, _p)
            log(f'{_MOSEK_LIC_FILE_ENV} not found but mosek file on another path:', p)
            break
# store for logging
MOSEK_LICENCE_FILE_PATH = os.getenv(_MOSEK_LIC_FILE_ENV)

safe_create_dir(DATA_DIR)


def get_safe_data_file(file_path: Path):
    """

    Parameters
    ----------
    file_path: :obj:`pathlib.Path`
    Data file path, if it is absolute will use as is, otherwise will look under data folder

    Returns
    -------
    normalized data file path

    """
    file_path = get_safe_path_obj(file_path)
    if file_path.is_absolute():
        return file_path
    return DATA_DIR / file_path
