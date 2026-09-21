"""Generic file download and tar.gz extraction helpers."""
import shutil
import tarfile
from urllib.request import urlopen


# Per-socket-operation timeout (not a cap on the whole transfer): a stalled
# connection raises instead of hanging the run forever.
DOWNLOAD_TIMEOUT_SECONDS = 120


def download(url, file_path):
    """Download `url` to `file_path`."""
    with urlopen(url, timeout=DOWNLOAD_TIMEOUT_SECONDS) as response, open(file_path, 'wb') as out_file:
        shutil.copyfileobj(response, out_file)


def extract_targz(targz_path, output_dir):
    """Extract a tar archive (compression auto-detected from its content,
    not its file name) into `output_dir`.

    Raises ValueError for anything that isn't a readable tar archive. This
    used to silently do nothing for a file not named *.tar.gz/*.tgz, which
    only surfaced later as a confusing missing-file error (no nodes.dmp, no
    UNITE files) far from the real cause.
    """
    try:
        archive = tarfile.open(str(targz_path), 'r:*')
    except tarfile.ReadError as e:
        raise ValueError(f'{targz_path} is not a readable tar archive ({e}).') from e
    with archive:
        # filter="data" rejects unsafe members (absolute paths, path traversal); required
        # explicitly pre-3.14 to avoid DeprecationWarning and to opt into the safer default.
        archive.extractall(path=output_dir, filter='data')
