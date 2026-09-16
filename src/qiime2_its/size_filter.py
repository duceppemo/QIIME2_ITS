"""Read-length filtering via BBDuk (BBTools)."""
import subprocess
from concurrent import futures
from pathlib import Path


def bbduk_version():
    """Returns `bbduk.sh --version`'s raw stderr text (BBTools writes its
    version banner there, not stdout), or None if bbduk.sh isn't installed
    -- it's optional, only needed for --min-len/--max-len. Used for the
    report's QA/provenance page (provenance.parse_bbduk_version() does the
    actual parsing).
    """
    try:
        result = subprocess.run(['bbduk.sh', '--version'], capture_output=True, text=True)
    except FileNotFoundError:
        return None
    return result.stderr


def size_select_se(fastq_in, out_folder, min_len, max_len, threads):
    out_path = Path(out_folder) / Path(fastq_in).name
    cmd = ['bbduk.sh', 'overwrite=t',
           f'in={fastq_in}', f'out={out_path}']
    if min_len > 0:
        cmd.append(f'minlength={min_len}')
    if max_len > 0:
        cmd.append(f'maxlength={max_len}')
    cmd.append(f'threads={threads}')
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)


def size_select_pe(fastq_r1, fastq_r2, out_folder, min_len, max_len, threads):
    out_r1 = Path(out_folder) / Path(fastq_r1).name
    out_r2 = Path(out_folder) / Path(fastq_r2).name
    cmd = ['bbduk.sh', 'overwrite=t',
           f'in={fastq_r1}', f'in2={fastq_r2}',
           f'out={out_r1}', f'out2={out_r2}']
    if min_len > 0:
        cmd.append(f'minlength={min_len}')
    if max_len > 0:
        cmd.append(f'maxlength={max_len}')
    cmd.append(f'threads={threads}')
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)


def size_select_se_parallel(fastq_list, out_folder, min_len, max_len, cpu, parallel):
    per_sample_cpu = max(1, int(cpu / parallel))
    with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
        list(executor.map(
            lambda fq: size_select_se(fq, out_folder, min_len, max_len, per_sample_cpu),
            fastq_list))


def size_select_pe_parallel(sample_dict, out_folder, min_len, max_len, cpu, parallel):
    per_sample_cpu = max(1, int(cpu / parallel))
    with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
        pairs = ((reads[0], reads[1]) for reads in sample_dict.values())
        list(executor.map(
            lambda p: size_select_pe(p[0], p[1], out_folder, min_len, max_len, per_sample_cpu),
            pairs))
