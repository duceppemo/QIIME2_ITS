"""Generic file download and tar.gz extraction helpers."""
import os
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
        if hasattr(tarfile, 'data_filter'):
            # filter="data" rejects unsafe members (absolute paths, path traversal, links
            # escaping output_dir, devices); explicit to opt into the safer default pre-3.14.
            archive.extractall(path=output_dir, filter='data')
        else:
            # Python 3.9-3.11 patch releases from before extraction filters existed
            # (< 3.9.17 / 3.10.12 / 3.11.4) raise TypeError on filter= -- apply the
            # essential part of the same policy by hand instead.
            root = os.path.realpath(output_dir)
            for member in archive.getmembers():
                target = os.path.realpath(os.path.join(root, member.name))
                if not (member.isfile() or member.isdir()) or os.path.commonpath([root, target]) != root:
                    raise ValueError(f'{targz_path} contains an unsafe member: {member.name}')
            archive.extractall(path=output_dir)
