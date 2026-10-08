#!/usr/bin/env python3
"""Generate docs/api.html (function reference) from the MATLAB sources.

Reads the COMMITTED version of every file (git HEAD by default), so the page
describes validated, committed code; pass --worktree to read the working tree.
For each function: signature, one-line summary, every argument with its size,
type, validators and default (from the `arguments` block), and the full help
text. Also the run scripts' headers and the legacy files in old/.

Run from anywhere:  python3 docs/tools/make_api.py
"""
import html
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]          # PulsarSimMatlab
WORKTREE = "--worktree" in sys.argv


# Hand-written descriptions of the legacy files (their headers are inconsistent).
LEGACY = {
    "applyDispersion.m": "In-memory dispersion by overlap-add. Superseded by <code>applyDispersionStream</code>; its kernel had "
                         "the wrong sign (<a href=\"lessons.html\">bug #1</a>).",
    "envelopeReconstruction.m": "Streamed IQ → envelope/power with optional decimation. Superseded by <code>detectPower</code>; "
                                "still uses the old bin timing and byte order.",
    "plotEnvelope.m": "Plots the output of <code>envelopeReconstruction</code>.",
    "process_large_baseband_file.m": "Early sketch: dynamic spectrum of an int8 complex baseband file in blocks.",
    "read_baseband_chunk.m": "Early sketch: read a segment of int8 complex samples.",
    "stream_pulsar_baseband_to_file.m": "Early sketch: stream 1 GS/s complex baseband pulsar data to disk.",
    "test.m": "Scratch check comparing raw and dispersed halves of <code>data/test.dat</code>.",
    "test_pipeline.m": "The first end-to-end testbed: filterbank-level simulation of several pulsars, dedispersion, folding, "
                       "TOAs and a position fix. Superseded by the voltage-level pipeline.",
}


def git(*args):
    return subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True,
                          text=True, check=True).stdout


def list_files(folder):
    if WORKTREE:
        return sorted(p.relative_to(ROOT).as_posix() for p in (ROOT / folder).glob("*.m"))
    out = git("ls-tree", "--name-only", "HEAD", folder + "/")
    return sorted(f for f in out.split() if f.endswith(".m"))


def read(rel):
    if WORKTREE:
        return (ROOT / rel).read_text()
    return git("show", "HEAD:./" + rel)


def help_block(src):
    """Text between the first %{ and %} (block comment), else leading % lines."""
    m = re.search(r"^%\{\s*$(.*?)^%\}\s*$", src, re.S | re.M)
    if m:
        return m.group(1).strip("\n")
    lines = []
    for ln in src.splitlines()[1:]:
        if ln.startswith("%"):
            lines.append(ln[1:])
        elif lines:
            break
    return "\n".join(lines)


def summary(src):
    for ln in src.splitlines()[:3]:
        m = re.match(r"%\s*[A-Z0-9_]+\s+(.*)", ln)
        if m:
            return m.group(1).strip()
    return ""


def arguments(src):
    m = re.search(r"^arguments\s*$(.*?)^end\s*$", src, re.S | re.M)
    if not m:
        return []
    rows = []
    for ln in m.group(1).splitlines():
        ln = ln.split("%")[0].rstrip()
        if not ln.strip():
            continue
        default = ""
        if "=" in ln:
            ln, default = ln.split("=", 1)
            default = default.strip()
        parts = ln.strip()
        name = parts.split()[0]
        rest = parts[len(name):].strip()
        size = ""
        if rest.startswith("("):
            size = rest[:rest.index(")") + 1]
            rest = rest[len(size):].strip()
        validators = ""
        if "{" in rest:
            validators = rest[rest.index("{"):].strip()
            rest = rest[:rest.index("{")].strip()
        rows.append((name, size, rest, validators, default))
    return rows


def esc(s):
    return html.escape(s, quote=False)


def code(s):
    return f"<code>{esc(s)}</code>" if s else ""


def anchor(name):
    return name.lower()


def function_section(rel):
    src = read(rel)
    first = src.splitlines()[0]
    name = Path(rel).stem
    sig = re.sub(r"^function\s+", "", first).strip()
    out = [f'<h3 id="{anchor(name)}"><code>{esc(name)}</code></h3>',
           f'<p>{esc(summary(src))}</p>',
           f'<div class="sig">{esc(sig)}</div>']
    rows = arguments(src)
    if rows:
        out.append('<div class="table-wrap"><table><tr><th>Argument</th><th>Size</th>'
                   '<th>Class</th><th>Validators</th><th>Default</th></tr>')
        for n, size, cls, val, d in rows:
            opt = n.startswith("opts.")
            label = f"'{n[5:]}'" if opt else n
            out.append(f"<tr><td><code>{esc(label)}</code></td><td>{esc(size)}</td><td>{esc(cls)}</td>"
                       f"<td>{code(val)}</td><td>{code(d) if d else ('' if opt else 'required')}</td></tr>")
        out.append("</table></div>")
    hb = help_block(src)
    if hb:
        out.append(f"<details><summary>Full help text</summary><pre><code>{esc(hb)}</code></pre></details>")
    return "\n".join(out)


def script_section(rel):
    src = read(rel)
    name = Path(rel).stem
    hb = help_block(src)
    first = summary(src) or (hb.strip().splitlines()[0] if hb.strip() else "")
    return "\n".join([f'<h3 id="{anchor(name)}"><code>{esc(Path(rel).name)}</code></h3>',
                      f"<p>{esc(first)}</p>",
                      f"<details><summary>Header</summary><pre><code>{esc(hb)}</code></pre></details>"])


def main():
    funcs = list_files("functions")
    scripts = [f for f in list_files(".") if "/" not in f]
    legacy = list_files("old")
    commit = git("rev-parse", "--short", "HEAD").strip()
    source = "working tree" if WORKTREE else f"commit <code>{commit}</code>"

    toc = " · ".join(f'<a href="#{anchor(Path(f).stem)}"><code>{Path(f).stem}</code></a>' for f in funcs)
    body = [f"<p class=\"meta\">Generated by <code>docs/tools/make_api.py</code> from {source}. "
            "Do not edit by hand: rerun the script after a code change.</p>",
            '<p class="lede">Every function in <code>functions/</code> with its signature, every argument and default '
            '(from the <code>arguments</code> block) and its full help text. The stage pages explain the <em>why</em>; '
            'this page is the exact interface. In signatures, <code>opts</code> stands for the name-value options listed '
            'in quotes.</p>',
            f'<div class="toc"><b>Functions</b><p>{toc}</p></div>',
            "<h2>Functions</h2>"]
    body += [function_section(f) for f in funcs]
    body.append("<h2>Scripts</h2>")
    body += [script_section(f) for f in scripts]
    body.append("<h2>Legacy (<code>old/</code>, unused)</h2>")
    body.append("<p>Superseded code kept for reference. Nothing calls it. <code>envelopeReconstruction</code> and "
                "<code>plotEnvelope</code> are on the open-items list (retire, or update to <code>binTime0</code> and "
                "<code>'ieee-le'</code>).</p>")
    body.append('<div class="table-wrap"><table><tr><th>File</th><th>Summary</th></tr>')
    for f in legacy:
        name = Path(f).name
        s = LEGACY.get(name)
        if s is None:
            src = read(f)
            s = esc(summary(src) or (help_block(src).strip().splitlines() or [""])[0])
        body.append(f"<tr><td><code>{esc(name)}</code></td><td>{s}</td></tr>")
    body.append("</table></div>")

    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Function Reference</title>
<link rel="stylesheet" href="style.css">
</head>
<body>
<main>

<h1>Function reference</h1>
{chr(10).join(body)}

</main>
<script src="nav.js"></script>
</body>
</html>
"""
    (ROOT / "docs" / "api.html").write_text(page)
    print(f"wrote docs/api.html: {len(funcs)} functions, {len(scripts)} scripts, {len(legacy)} legacy files")


if __name__ == "__main__":
    main()
