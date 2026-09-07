#!/usr/bin/env python3
"""Pack a golden fixture directory into one tracked text file, or unpack it.

    ./gen/pack_golden.py pack   tests/golden/matrix
    ./gen/pack_golden.py unpack tests/golden/matrix.pack

The fixtures were 1914 loose files, 95% of everything this repository tracks,
and git pays for every one of them on every checkout and every `status`. They
are 0.5 MB of text between them, so the file count was the whole cost and none
of the value.

WHY A TEXT CONTAINER AND NOT A TAR

A tarball would collapse them just as well and would be the obvious choice, but
it would cost the property these fixtures exist for. `tests/matrix.rs` puts it
plainly: "a fixture that changes is the port changing its answer, which is the
thing this repository exists to prevent." That only works if a reviewer can see
the change, and a diff of a compressed archive is one unreadable blob. A text
container diffs line by line, so a fixture that moves still shows up as moved
lines in review, and git packs the result to about what the tarball would be.

THE FORMAT

    # mafft-rs golden pack v1
    <relative path> <byte count>
    <exactly that many bytes>
    <relative path> <byte count>
    ...

Length-prefixed rather than delimited, so no fixture content can imitate a
record boundary -- these files hold arbitrary alignment output, and any
sentinel line would eventually appear inside one of them.
"""

import os
import sys
from pathlib import Path

HEADER = b"# mafft-rs golden pack v1\n"


def pack(src: Path) -> Path:
    """Write `<src>.pack` from every file under `src`, sorted for a stable diff."""
    out = src.with_suffix(".pack")
    files = sorted(p for p in src.rglob("*") if p.is_file())
    with out.open("wb") as fh:
        fh.write(HEADER)
        for p in files:
            rel = p.relative_to(src).as_posix()
            data = p.read_bytes()
            fh.write(f"{rel} {len(data)}\n".encode())
            fh.write(data)
            if data and not data.endswith(b"\n"):
                # Keeps the next header on its own line without changing the
                # payload: the byte count, not this newline, is what is read.
                fh.write(b"\n")
    return out


def unpack(pack_file: Path, dst: Path) -> int:
    """Recreate the directory `pack_file` was made from. Returns the file count."""
    raw = pack_file.read_bytes()
    if not raw.startswith(HEADER):
        raise SystemExit(f"{pack_file}: not a golden pack")
    i = len(HEADER)
    n = 0
    while i < len(raw):
        nl = raw.index(b"\n", i)
        rel, count = raw[i:nl].decode().rsplit(" ", 1)
        i = nl + 1
        data = raw[i : i + int(count)]
        i += int(count)
        if i < len(raw) and raw[i : i + 1] == b"\n":
            i += 1
        target = dst / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
        n += 1
    return n


def main() -> None:
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    mode, arg = sys.argv[1], Path(sys.argv[2])
    if mode == "pack":
        out = pack(arg)
        print(f"{out}: {os.path.getsize(out)} bytes")
    elif mode == "unpack":
        dst = Path(sys.argv[3]) if len(sys.argv) > 3 else arg.with_suffix("")
        print(f"{dst}: {unpack(arg, dst)} files")
    else:
        raise SystemExit(__doc__)


main()
