"""Generic file download and tar.gz extraction helpers."""
import shutil
import tarfile
from urllib.request import urlopen


def download(url, file_path):
    """Download `url` to `file_path`."""
    with urlopen(url) as response, open(file_path, 'wb') as out_file:
        shutil.copyfileobj(response, out_file)


def extract_targz(targz_path, output_dir):
    """Extract a .tar.gz/.tgz archive into `output_dir`."""
    targz_path = str(targz_path)
    if not targz_path.endswith(('.tar.gz', '.tgz')):
        return
    with tarfile.open(targz_path, 'r:gz') as archive:
        # filter="data" rejects unsafe members (absolute paths, path traversal); required
        # explicitly pre-3.14 to avoid DeprecationWarning and to opt into the safer default.
        archive.extractall(path=output_dir, filter='data')
