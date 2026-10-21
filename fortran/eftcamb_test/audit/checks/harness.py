"""EFTCAMB audit harness: run any ini on any built tree, save outputs to npz, compare two runs.

Every run happens in its own subprocess (crash isolation, own OMP_NUM_THREADS, cwd = output directory), importing
camb from the selected tree only.  Trees are the frozen builds made by setup_trees.py under ``$AUDIT_ROOT/trees``.

Command line:
    python harness.py run   --tree A0 --ini MODEL.ini --out RUN.npz [--set key=value ...] [--evol-k 0.01,0.1]
    python harness.py batch --tree A0 --inis 'parameters/*.ini' --outdir DIR [--set key=value ...] [--jobs 8]
    python harness.py compare A.npz B.npz [--json]

Library use (from other audit scripts):
    from harness import run_model, run_batch, compare, load, model_inis, TREES, slot, build_tree

Resource limits (shared by every process on the machine that uses this module):
    CPU work is gated by a machine-wide pool of AUDIT_MAX_SLOTS (default 8) thread slots held through file locks under
    $AUDIT_ROOT/locks; run_model takes ``threads`` slots, ``with slot(n):`` guards any other CAMB work, and
    build_tree serialises tree rebuilds (one at a time).  Workers run at nice 10 and exit if their parent
    process dies.  Locks are released automatically when a holder exits.
"""

from __future__ import annotations

import argparse
import contextlib
import fcntl
import glob
import hashlib
import json
import os
import shutil
import signal
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
PARAMS_DIR = REPO / "fortran" / "eftcamb_test" / "parameters"
DEFAULT_ROOT = (
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
    "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
)
AUDIT_ROOT = Path(os.environ.get("AUDIT_ROOT", DEFAULT_ROOT))
TREES = {name: AUDIT_ROOT / "trees" / name for name in ("A0", "AD", "ADD", "B0", "B1")}

# default grids saved with every run
A_GRID = np.logspace(-8, 0, 401)
Z_BG = np.concatenate([np.linspace(0, 10, 201), np.logspace(1.05, 4, 60)])
CL_KEYS = ("total", "unlensed_scalar", "unlensed_total", "lensed_scalar", "tensor", "lens_potential")

MAX_SLOTS = int(os.environ.get("AUDIT_MAX_SLOTS", "8"))
LOCK_DIR = AUDIT_ROOT / "locks"


@contextlib.contextmanager
def slot(n: int = 1, poll: float = 0.5):
    """Hold ``n`` of the machine-wide CPU slots (blocking) for the duration of the block."""
    n = max(1, min(int(n), MAX_SLOTS))
    LOCK_DIR.mkdir(parents=True, exist_ok=True)
    held: list = []
    try:
        while True:
            for i in range(MAX_SLOTS):
                if len(held) == n:
                    break
                f = open(LOCK_DIR / f"cpu_{i}.lock", "w")
                try:
                    fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    held.append(f)
                except BlockingIOError:
                    f.close()
            if len(held) == n:
                break
            # release partial holdings so that larger requests cannot deadlock each other
            for f in held:
                f.close()
            held = []
            time.sleep(poll)
        yield
    finally:
        for f in held:
            f.close()


@contextlib.contextmanager
def _build_lock():
    LOCK_DIR.mkdir(parents=True, exist_ok=True)
    with open(LOCK_DIR / "build.lock", "w") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        yield


def _nice() -> None:
    os.nice(10)


def build_tree(src: str | Path, dest: str | Path, fflags: str | None = None, make_args: list[str] | None = None) -> dict:
    """Copy a tree (if ``dest`` does not exist) and (re)build it: one build at a time machine-wide, nice 10.

    fflags=None builds with release flags via setup.py make; otherwise make python FFLAGS=<fflags> in fortran/.
    Edit files in ``dest`` between calls and call again to rebuild.
    """
    src, dest = tree_path(src), Path(dest)
    if not dest.exists():
        shutil.copytree(src, dest, symlinks=True)
    # build artefacts live in the scratch area, keyed by the tree path, never inside the tree
    key = f"{dest.name}_{hashlib.md5(str(dest.resolve()).encode()).hexdigest()[:8]}"
    env = dict(os.environ, CAMB_BUILD_DIR=str(AUDIT_ROOT / "build" / key))
    if fflags is None and not make_args:
        cmd, cwd = [sys.executable, "setup.py", "make"], dest
    else:
        cmd, cwd = ["make", "python", "PYCAMB_OUTPUT_DIR=../camb/", *(make_args or [])], dest / "fortran"
        if fflags is not None:
            cmd.append(f"FFLAGS={fflags}")
    lib = dest / "camb" / "camblib.so"
    mtime_before = lib.stat().st_mtime if lib.exists() else 0.0
    log = AUDIT_ROOT / "logs" / f"build_{key}.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with _build_lock(), slot(4), open(log, "w") as f:
        rc = subprocess.run(cmd, cwd=cwd, env=env, stdout=f, stderr=subprocess.STDOUT, preexec_fn=_nice).returncode
    # setup.py make can exit 0 on a failed compile: also check the log (fresh_library tells whether it was rewritten)
    text = log.read_text(errors="replace")
    fresh = lib.exists() and lib.stat().st_mtime > mtime_before
    tail = "\n".join(text.splitlines()[-50:])
    ok = rc == 0 and lib.exists() and "Compilation failed" not in text and "Error " not in tail
    return {"ok": ok, "fresh_library": fresh, "returncode": rc, "tree": str(dest), "log": str(log), "elapsed": round(time.time() - t0, 1)}


def model_inis(pattern: str = "*.ini") -> list[Path]:
    """Model inis of the EFTCAMB test suite (excludes the base parameter files)."""
    return sorted(p for p in PARAMS_DIR.glob(pattern) if "base_params" not in p.name)


def tree_path(tree: str | Path) -> Path:
    return TREES[tree] if str(tree) in TREES else Path(tree)


def _wrapper_ini(ini: Path, outdir: Path, name: str, overrides: dict[str, str]) -> Path:
    lines = [
        f"DEFAULT({ini.resolve()})",
        f"output_root = {outdir / name}",
        "EFTCAMB_write_background = F",
    ]
    lines += [f"{k} = {v}" for k, v in overrides.items()]
    path = outdir / f"{name}.wrapper.ini"
    path.write_text("\n".join(lines) + "\n")
    return path


def run_model(
    tree: str | Path,
    ini: str | Path,
    out: str | Path,
    overrides: dict[str, str] | None = None,
    evol_k: list[float] | None = None,
    threads: int = 1,
    timeout: float = 1800,
    extra_env: dict[str, str] | None = None,
    python_set: dict[str, object] | None = None,
) -> dict:
    """Run one model in a subprocess. Returns a status dict (also written next to ``out`` as .json).

    overrides: ini key -> value, appended after DEFAULT(model ini) so they take precedence.
    python_set: attribute path -> value set on CAMBparams after read_ini (e.g. {"EFTCAMB.EFTCAMB_GR_threshold": 1e-9}).
    """
    out = Path(out).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    ini = Path(ini)
    name = out.stem
    wrapper = _wrapper_ini(ini, out.parent, name, overrides or {})
    cfg = {
        "tree": str(tree_path(tree)),
        "ini": str(wrapper),
        "model_ini": str(ini.resolve()),
        "out": str(out),
        "evol_k": evol_k or [],
        "python_set": python_set or {},
    }
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), PYTHONHASHSEED="0")
    env.pop("PYTHONPATH", None)
    env.update(extra_env or {})
    log = out.with_suffix(".log")
    t_wait = time.time()
    try:
        with slot(threads), open(log, "w") as f:
            t0 = time.time()
            proc = subprocess.run(
                [sys.executable, str(HERE / "harness.py"), "_worker", json.dumps(cfg)],
                cwd=out.parent,
                env=env,
                stdout=f,
                stderr=subprocess.STDOUT,
                timeout=timeout,
                preexec_fn=_nice,
            )
        rc = proc.returncode
        status = "ok" if rc == 0 else (f"signal:{signal.Signals(-rc).name}" if rc < 0 else f"error:rc={rc}")
    except subprocess.TimeoutExpired:
        rc, status = None, "timeout"
        t0 = t_wait
    info = {
        "tree": str(tree),
        "model_ini": str(ini),
        "overrides": overrides or {},
        "python_set": python_set or {},
        "threads": threads,
        "status": status,
        "returncode": rc,
        "elapsed": round(time.time() - t0, 2),
        "out": str(out),
        "log": str(log),
    }
    if status != "ok":
        info["log_tail"] = log.read_text(errors="replace").splitlines()[-40:]
    out.with_suffix(".json").write_text(json.dumps(info, indent=2))
    return info


def run_batch(
    tree: str | Path, inis: list[Path], outdir: str | Path, jobs: int = 6, suffix: str = "", **kw
) -> list[dict]:
    outdir = Path(outdir)

    def one(ini: Path) -> dict:
        return run_model(tree, ini, outdir / f"{ini.stem}{suffix}.npz", **kw)

    with ThreadPoolExecutor(jobs) as ex:
        return list(ex.map(one, inis))


def load(path: str | Path) -> dict:
    with np.load(path, allow_pickle=False) as f:
        d = {k: f[k] for k in f.files}
    d["meta"] = json.loads(str(d["meta"])) if "meta" in d else {}
    return d


# --------------------------------------------------------------------------------------------------
# comparison


def _rel(a: np.ndarray, b: np.ndarray, floor: float = 0.0) -> float:
    den = np.maximum(np.abs(b), floor)
    m = den > 0
    return float(np.max(np.abs(a[m] - b[m]) / den[m])) if m.any() else 0.0


def compare_cls(a: np.ndarray, b: np.ndarray, lmin: int = 2, lmax: int = 2500) -> dict:
    """Max fractional differences for a (lmax+1, 4) TT, EE, BB, TE array; TE normalised by sqrt(TT EE)."""
    n = min(len(a), len(b), lmax + 1)
    a, b = a[lmin:n], b[lmin:n]
    res = {}
    for i, lab in enumerate(("TT", "EE", "BB")):
        if np.any(b[:, i] != 0):
            res[lab] = _rel(a[:, i], b[:, i])
    if a.shape[1] > 3:
        norm = np.sqrt(np.abs(b[:, 0] * b[:, 1]))
        m = norm > 0
        res["TE"] = float(np.max(np.abs(a[m, 3] - b[m, 3]) / norm[m])) if m.any() else 0.0
    return res


def compare(a: dict | str | Path, b: dict | str | Path, lmax: int = 2500) -> dict:
    """Compare two harness outputs; returns max fractional differences per observable and bit-identity flags."""
    A = load(a) if not isinstance(a, dict) else a
    B = load(b) if not isinstance(b, dict) else b
    res: dict = {"bit_identical": True, "cl": {}, "missing": []}
    for key in CL_KEYS:
        ka = f"cl_{key}"
        if ka in A and ka in B:
            if key == "lens_potential":
                n = min(len(A[ka]), len(B[ka]), lmax + 1)
                res["cl"][key] = {"PP": _rel(A[ka][2:n, 0], B[ka][2:n, 0])}
            else:
                res["cl"][key] = compare_cls(A[ka], B[ka], lmax=lmax)
            n = min(len(A[ka]), len(B[ka]))
            res["bit_identical"] &= bool(np.array_equal(A[ka][:n], B[ka][:n]))
        elif ka in A or ka in B:
            res["missing"].append(ka)
    if "pk" in A and "pk" in B and A["pk"].shape == B["pk"].shape and np.allclose(A["pk_k"], B["pk_k"], rtol=1e-12):
        res["pk"] = [_rel(A["pk"][i], B["pk"][i]) for i in range(A["pk"].shape[0])]
        res["bit_identical"] &= bool(np.array_equal(A["pk"], B["pk"]))
    elif "pk" in A and "pk" in B:
        kmin = max(A["pk_k"][0], B["pk_k"][0])
        kmax = min(A["pk_k"][-1], B["pk_k"][-1])
        k = np.logspace(np.log10(kmin), np.log10(kmax), 400)
        nz = min(A["pk"].shape[0], B["pk"].shape[0])
        res["pk"] = [
            _rel(
                np.exp(np.interp(np.log(k), np.log(A["pk_k"]), np.log(A["pk"][i]))),
                np.exp(np.interp(np.log(k), np.log(B["pk_k"]), np.log(B["pk"][i]))),
            )
            for i in range(nz)
        ]
        res["bit_identical"] = False
    for key in ("sigma8", "bg_H", "bg_tau"):
        if key in A and key in B and A[key].shape == B[key].shape:
            res[key] = _rel(A[key], B[key])
            res["bit_identical"] &= bool(np.array_equal(A[key], B[key]))
    eft_a = {k[4:] for k in A if k.startswith("eft_")}
    eft_b = {k[4:] for k in B if k.startswith("eft_")}
    if eft_a & eft_b:
        res["eft"] = {}
        for f in sorted(eft_a & eft_b):
            x, y = A[f"eft_{f}"], B[f"eft_{f}"]
            if x.shape == y.shape and x.dtype.kind == "f":
                scale = np.max(np.abs(y)) if np.any(y) else 1.0
                res["eft"][f] = float(np.max(np.abs(x - y)) / scale) if scale else 0.0
    res["worst_cl"] = max([v for d in res["cl"].values() for v in d.values()] or [0.0])
    res["worst_pk"] = max(res.get("pk", [0.0]) or [0.0])
    return res


# --------------------------------------------------------------------------------------------------
# worker (runs inside the subprocess)


def parse_ini(path: str | Path) -> dict[str, str]:
    """Flatten a CAMB ini file: keys in a file override those from its DEFAULT(...) files; INCLUDE(...) is inline."""
    path = Path(path)
    keys: dict[str, str] = {}
    defaults: list[dict[str, str]] = []
    for raw in path.read_text().splitlines():
        line = raw.split("#", 1)[0].strip() if not raw.lstrip().startswith("#") else ""
        if not line:
            continue
        up = line.upper()
        if up.startswith("DEFAULT(") or up.startswith("INCLUDE("):
            inc = Path(line[line.index("(") + 1 : line.rindex(")")].strip())
            inc = inc if inc.is_absolute() else path.parent / inc
            (defaults.append if up.startswith("DEFAULT(") else keys.update)(parse_ini(inc))
        elif "=" in line:
            k, v = line.split("=", 1)
            keys[k.strip()] = v.strip()
    merged: dict[str, str] = {}
    for d in defaults:
        merged.update(d)
    merged.update(keys)
    return merged


def _set_attr(obj, path: str, value) -> None:
    parts = path.split(".")
    for p in parts[:-1]:
        obj = getattr(obj, p)
    setattr(obj, parts[-1], value)


def _exit_with_parent() -> None:
    parent = os.getppid()

    def watch() -> None:
        while True:
            time.sleep(2)
            if os.getppid() != parent:
                os._exit(3)

    threading.Thread(target=watch, daemon=True).start()


def _worker(cfg: dict) -> None:
    _exit_with_parent()
    sys.path.insert(0, cfg["tree"])
    import camb

    assert Path(camb.__file__).resolve().is_relative_to(Path(cfg["tree"]).resolve()), camb.__file__
    out: dict[str, np.ndarray] = {}
    meta: dict = {"camb_file": camb.__file__, "camb_version": getattr(camb, "__version__", "?"), "errors": {}}
    t0 = time.time()
    cp = camb.read_ini(cfg["ini"], no_validate=True)
    for path, value in cfg["python_set"].items():
        _set_attr(cp, path, value)
    has_eft = hasattr(cp, "EFTCAMB")
    ini_keys = parse_ini(cfg["ini"])
    if has_eft and cp.EFTCAMB.EFTflag == 0 and int(ini_keys.get("EFTflag", "0")) != 0:
        # pre-rebase trees (B1) do not read EFTCAMB keys in read_ini: initialise EFTCAMB from the flattened ini
        cp.EFTCAMB.initialize_parameters(cp, ini_keys, print_header=False)
        meta["eft_init_fallback"] = True
    # The character(len=200) getters (model_name, read_parameters, param_names) are called from Python without the
    # hidden gfortran length argument, which -fcheck=bounds rejects; skip them on the runtime-check trees.
    string_calls = cfg.get("string_calls", Path(cfg["tree"]).name not in ("AD", "ADD"))
    meta["string_calls"] = string_calls
    if has_eft and not string_calls:
        meta["EFTflag"] = int(cp.EFTCAMB.EFTflag)
    if has_eft and string_calls:
        try:
            meta["eft_params"] = {k: (v if isinstance(v, (int, float, bool)) else str(v))
                                  for k, v in cp.EFTCAMB.read_parameters().items()}
            meta["eft_model_name"] = cp.EFTCAMB.model_name()
            meta["EFTflag"] = int(cp.EFTCAMB.EFTflag)
        except Exception as e:  # noqa: BLE001
            meta["errors"]["eft_meta"] = repr(e)
    t1 = time.time()
    results = camb.get_results(cp)
    meta["time_get_results"] = time.time() - t1
    meta["time_read"] = t1 - t0
    try:
        meta["derived"] = {k: float(v) for k, v in results.get_derived_params().items()}
    except Exception as e:  # noqa: BLE001
        meta["errors"]["derived"] = repr(e)
    if cp.WantCls:
        try:
            cls = results.get_cmb_power_spectra(CMB_unit="muK", raw_cl=True)
            for key in CL_KEYS:
                if key in cls:
                    out[f"cl_{key}"] = np.asarray(cls[key])
        except Exception as e:  # noqa: BLE001
            meta["errors"]["cls"] = repr(e)
    if cp.WantTransfer:
        try:
            k, z, pk = results.get_linear_matter_power_spectrum(hubble_units=False, k_hunit=False)
            out["pk_k"], out["pk_z"], out["pk"] = np.asarray(k), np.asarray(z), np.asarray(pk)
            out["sigma8"] = np.atleast_1d(results.get_sigma8())
        except Exception as e:  # noqa: BLE001
            meta["errors"]["pk"] = repr(e)
    try:
        out["bg_z"] = Z_BG
        out["bg_H"] = np.asarray(results.hubble_parameter(Z_BG))
        out["bg_tau"] = np.asarray(results.conformal_time(Z_BG))
        out["bg_DA"] = np.asarray(results.angular_diameter_distance(Z_BG))
    except Exception as e:  # noqa: BLE001
        meta["errors"]["background"] = repr(e)
    if has_eft and cp.EFTCAMB.EFTflag != 0:
        try:
            fields, arr = cp.EFTCAMB.get_eft_functions(results, A_GRID)
            out["eft_a_grid"] = A_GRID
            for f in arr.dtype.names:
                out[f"eft_{f}"] = np.ascontiguousarray(arr[f])
            meta["eft_fields"] = list(arr.dtype.names)
        except Exception as e:  # noqa: BLE001
            meta["errors"]["eft_functions"] = repr(e)
        if cfg["evol_k"]:
            try:
                a_ev = np.logspace(-6, 0, 300)
                fields, arr = cp.EFTCAMB.get_scale_evolution(results, cfg["evol_k"], a_ev)
                out["evol_k"], out["evol_a"] = np.asarray(cfg["evol_k"]), a_ev
                for f in arr.dtype.names:
                    out[f"evol_{f}"] = np.ascontiguousarray(arr[f])
            except Exception as e:  # noqa: BLE001
                meta["errors"]["evolution"] = repr(e)
    meta["time_total"] = time.time() - t0
    out["meta"] = np.array(json.dumps(meta))
    np.savez_compressed(cfg["out"], **out)
    print(json.dumps({"ok": True, "errors": meta["errors"]}))


def _parse_set(items: list[str]) -> dict[str, str]:
    return dict(s.split("=", 1) for s in items)


def main() -> None:
    if len(sys.argv) > 2 and sys.argv[1] == "_worker":
        _worker(json.loads(sys.argv[2]))
        return
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    bd = sub.add_parser("build", help="build a tree in place (machine-wide build lock, nice 10)")
    bd.add_argument("--tree", required=True)
    bd.add_argument("--fflags", default=None)
    r = sub.add_parser("run")
    r.add_argument("--tree", required=True)
    r.add_argument("--ini", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--set", nargs="*", default=[])
    r.add_argument("--evol-k", default="")
    r.add_argument("--threads", type=int, default=2)
    b = sub.add_parser("batch")
    b.add_argument("--tree", required=True)
    b.add_argument("--inis", default="*.ini")
    b.add_argument("--outdir", required=True)
    b.add_argument("--set", nargs="*", default=[])
    b.add_argument("--jobs", type=int, default=6)
    b.add_argument("--threads", type=int, default=2)
    c = sub.add_parser("compare")
    c.add_argument("a")
    c.add_argument("b")
    c.add_argument("--json", action="store_true")
    a = p.parse_args()
    if a.cmd == "build":
        res = build_tree(a.tree, a.tree, fflags=a.fflags)
        print(json.dumps(res, indent=2))
        sys.exit(0 if res["ok"] else 1)
    if a.cmd == "run":
        evk = [float(x) for x in a.evol_k.split(",") if x]
        print(json.dumps(run_model(a.tree, a.ini, a.out, _parse_set(a.set), evk, a.threads), indent=2))
    elif a.cmd == "batch":
        inis = [Path(x) for x in sorted(glob.glob(a.inis))] if os.sep in a.inis else model_inis(a.inis)
        res = run_batch(a.tree, inis, a.outdir, a.jobs, overrides=_parse_set(a.set), threads=a.threads)
        for r_ in res:
            print(f"{Path(r_['model_ini']).stem:40s} {r_['status']:20s} {r_['elapsed']:8.1f}s")
    else:
        res = compare(a.a, a.b)
        print(json.dumps(res, indent=2) if a.json else f"worst_cl={res['worst_cl']:.3e} worst_pk={res['worst_pk']:.3e} "
              f"bit_identical={res['bit_identical']}")


if __name__ == "__main__":
    main()
