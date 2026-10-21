"""Shared file-format checks; no tool execution or workflow orchestration."""
import csv
import gzip
import hashlib
import json
import re
from collections import OrderedDict
from pathlib import Path


# Match init_mask_table() in Cactus v3.3.0's
# preprocessor/cactus_sanitizeFastaHeaders.c; uppercase derived inputs as before.
CACTUS_MASK_BASES = b'UWSMKRYBDHV'
CACTUS_MASK_TABLE = bytes.maketrans(CACTUS_MASK_BASES, b'N' * len(CACTUS_MASK_BASES))
CACTUS_BASE_POLICY = 'uppercase; UWSMKRYBDHV to N (Cactus 3.3.0); coordinates unchanged'


def cactus_bases(bases, context):
    """Return normalized bytes and masked-base count; reject non-Cactus symbols."""
    if not bases:
        raise ValueError(f'Empty sequence: {context}')
    bases = bases.upper()
    ambiguous = bases.translate(None, b'ACGTN')
    invalid = ambiguous.translate(None, CACTUS_MASK_BASES)
    if invalid:
        raise ValueError(f'Unsupported sequence symbols {bytes(sorted(set(invalid)))!r}: {context}')
    return bases.translate(CACTUS_MASK_TABLE), len(ambiguous)


def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def rows(path):
    with open_text(path) as handle:
        for line in handle:
            if line.strip():
                yield line.rstrip('\r\n').split('\t')


def table(path):
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_table(path, fields, data):
    with open(path, 'w') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(data)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def new_output_dir(path):
    path = Path(path)
    if path.exists() and any(path.iterdir()):
        raise ValueError(f'Refusing to append to a nonempty output directory: {path}. Resume through Snakemake.')
    path.mkdir(parents=True, exist_ok=True)
    return path


def identity(path, checksum=False):
    path = Path(path)
    stat = path.stat()
    if not path.is_file():
        raise ValueError(f'Not a file: {path}')
    result = dict(path=str(path), size=stat.st_size, mtime_ns=stat.st_mtime_ns)
    if checksum:
        digest = hashlib.sha256()
        with open(path, 'rb') as handle:
            for block in iter(lambda: handle.read(1048576), b''):
                digest.update(block)
        result['sha256'] = digest.hexdigest()
    return result


def unchanged(record):
    current = identity(record['path'], 'sha256' in record)
    if current != record:
        raise ValueError(f"Frozen input changed: {record['path']}; use a new build ID")


def check_paf(row):
    if len(row) < 13 or row[4] not in ('+', '-'):
        raise ValueError('Malformed node-level PAF record')
    qlen, qs, qe, tlen, ts, te = map(int, [row[1], row[2], row[3], row[6], row[7], row[8]])
    if not (0 <= qs < qe <= qlen and 0 <= ts < te <= tlen):
        raise ValueError(f'Invalid PAF coordinates: {row[:9]}')
    matches, block, mapq = map(int, row[9:12])
    if not 0 <= matches <= block or block <= 0 or not 0 <= mapq <= 255:
        raise ValueError('Invalid PAF match count, block length or mapping quality')
    cigars = [s[5:] for s in row[12:] if s.startswith('cg:Z:')]
    if len(cigars) != 1 or not re.fullmatch(r'(?:[1-9][0-9]*[MID=X])+', cigars[0]):
        raise ValueError(f'Expected one unclipped MID=X CIGAR: {row[0]}')
    ops = [(int(n), op) for n, op in re.findall(r'(\d+)([MID=X])', cigars[0])]
    if sum(n for n, op in ops if op != 'D') != qe - qs or sum(n for n, op in ops if op != 'I') != te - ts:
        raise ValueError(f'CIGAR consumption differs from PAF coordinates: {row[0]}')
    return qs, qe


class Writers:
    """Bound file descriptors even when the seed has many isolated components."""
    def __init__(self, limit=32):
        self.handles = OrderedDict()
        self.limit = limit

    def get(self, path):
        path = Path(path)
        if path not in self.handles:
            if len(self.handles) >= self.limit:
                self.handles.popitem(last=False)[1].close()
            path.parent.mkdir(parents=True, exist_ok=True)
            self.handles[path] = open(path, 'a')
        self.handles.move_to_end(path)
        return self.handles[path]

    def close(self):
        for handle in self.handles.values():
            handle.close()
        self.handles.clear()
