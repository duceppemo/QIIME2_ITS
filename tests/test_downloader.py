import tarfile

from qiime2_its import downloader


def test_download_writes_response_body_to_file(tmp_path, mocker):
    fake_response = mocker.MagicMock()
    fake_response.__enter__ = mocker.Mock(return_value=fake_response)
    fake_response.__exit__ = mocker.Mock(return_value=False)
    fake_response.read = mocker.Mock(side_effect=[b'hello ', b'world', b''])

    mocker.patch('qiime2_its.downloader.urlopen', return_value=fake_response)

    out_file = tmp_path / 'out.txt'
    downloader.download('http://example.org/file.txt', out_file)

    assert out_file.read_bytes() == b'hello world'


def test_extract_targz_extracts_members(tmp_path):
    src_dir = tmp_path / 'src'
    src_dir.mkdir()
    (src_dir / 'content.txt').write_text('payload')

    archive_path = tmp_path / 'archive.tar.gz'
    with tarfile.open(archive_path, 'w:gz') as tar:
        tar.add(src_dir / 'content.txt', arcname='content.txt')

    dest_dir = tmp_path / 'dest'
    dest_dir.mkdir()
    downloader.extract_targz(archive_path, dest_dir)

    assert (dest_dir / 'content.txt').read_text() == 'payload'


def test_extract_targz_ignores_non_archive_extension(tmp_path):
    not_an_archive = tmp_path / 'plain.txt'
    not_an_archive.write_text('not a tarball')
    dest_dir = tmp_path / 'dest'
    dest_dir.mkdir()

    downloader.extract_targz(not_an_archive, dest_dir)  # should not raise

    assert list(dest_dir.iterdir()) == []
