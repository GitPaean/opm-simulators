"""Minimal reader for ECLIPSE binary (Fortran-unformatted) output files."""
import numpy as np
import struct

_TYPE = {
    b'INTE': ('i4', 4, 1000),
    b'REAL': ('f4', 4, 1000),
    b'DOUB': ('f8', 8, 1000),
    b'LOGI': ('i4', 4, 1000),
    b'CHAR': ('S8', 8, 105),
    b'MESS': (None, 0, 1000),
}
for n in range(1, 100):
    _TYPE[b'C%03d' % n] = ('S%d' % n, n, 105)


def _rec(f):
    head = f.read(4)
    if len(head) < 4:
        return None
    (n,) = struct.unpack('>i', head)
    data = f.read(n)
    (n2,) = struct.unpack('>i', f.read(4))
    assert n == n2, 'record marker mismatch'
    return data


def read(path):
    """Yield (keyword, numpy array) in file order."""
    out = []
    with open(path, 'rb') as f:
        while True:
            hdr = _rec(f)
            if hdr is None:
                break
            kw = hdr[0:8].decode('ascii').strip()
            (num,) = struct.unpack('>i', hdr[8:12])
            typ = hdr[12:16]
            dtype, size, blk = _TYPE[typ]
            if dtype is None or num == 0:
                out.append((kw, np.array([])))
                continue
            chunks = []
            left = num
            while left > 0:
                take = min(left, blk)
                raw = _rec(f)
                if dtype.startswith('S'):
                    chunks.append(np.frombuffer(raw, dtype='S%d' % size, count=take))
                else:
                    chunks.append(np.frombuffer(raw, dtype='>' + dtype, count=take))
                left -= take
            out.append((kw, np.concatenate(chunks)))
    return out


def steps(path):
    """Split a unified restart file into {seqnum: {kw: array}}."""
    res, cur, seq = {}, None, None
    for kw, arr in read(path):
        if kw == 'SEQNUM':
            seq = int(arr[0])
            cur = {}
            res[seq] = cur
        elif cur is not None:
            cur.setdefault(kw, []).append(arr)
    return {s: {k: (v[0] if len(v) == 1 else np.concatenate(v)) for k, v in d.items()}
            for s, d in res.items()}


def summary(smspec, unsmry):
    """Return (dict name->series, times) from a SMSPEC/UNSMRY pair."""
    spec = dict()
    for kw, arr in read(smspec):
        spec.setdefault(kw, arr)
    names = [b.decode().strip() for b in spec['KEYWORDS']]
    wgn = [b.decode().strip() for b in spec.get('WGNAMES', spec.get('NAMES', []))]
    nums = spec.get('NUMS', np.zeros(len(names), dtype=int))
    labels = []
    for i, n in enumerate(names):
        w = wgn[i] if i < len(wgn) else ''
        if w and not w.startswith(':+:'):
            labels.append('%s:%s' % (n, w))
        elif int(nums[i]) > 0 and n.startswith(('B', 'R')):
            labels.append('%s:%d' % (n, int(nums[i])))
        else:
            labels.append(n)
    rows = [arr for kw, arr in read(unsmry) if kw == 'PARAMS']
    data = np.array(rows) if rows else np.zeros((0, len(labels)))
    return labels, data
