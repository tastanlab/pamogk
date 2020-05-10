from lib.sutils import *

ROOT_DIR = Path(__file__).resolve().parent
LIB_DIR = ROOT_DIR / 'lib'
DATA_DIR = ROOT_DIR / 'data'

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
