"""Adversarial verification of finding eft-par-cache-realloc-crash.

Each scenario runs in its own subprocess (a Fortran runtime error kills the interpreter).
Usage:
  python verify_par_cache_realloc.py scenarios --tree A0 --ini 4_Des_F_R_LCDM_1
  python verify_par_cache_realloc.py allmodels --tree A0 --jobs 4     (get_BAO(z) default on all test inis)
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(
    "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit"
)
PARAMS = Path(__file__).resolve().parents[3] / "parameters"
OUT = ROOT / "work" / "verify" / "par_cache_realloc"
Z = [0.5, 1.0]

CHILD = r"""
import json, sys
sys.path.insert(0, sys.argv[1])
sys.path.insert(0, sys.argv[5])
import camb, harness
tree, ini, mode = sys.argv[1], sys.argv[2], sys.argv[3]
Z = json.loads(sys.argv[4])
cp = camb.read_ini(ini, no_validate=True)
keys = harness.parse_ini(ini)
if hasattr(cp, "EFTCAMB") and cp.EFTCAMB.EFTflag == 0 and int(keys.get("EFTflag", "0")) != 0:
    cp.EFTCAMB.initialize_parameters(cp, keys, print_header=False)
flag = int(cp.EFTCAMB.EFTflag) if hasattr(cp, "EFTCAMB") else -1
print("EFTFLAG", flag, flush=True)
def H(r):
    return [float(x) for x in r.hubble_parameter(Z)]
if mode == "full_get_BAO":
    res = camb.get_results(cp)
else:
    res = camb.get_background(cp)
if mode == "fresh":
    out = H(res)
elif mode in ("get_BAO_default", "full_get_BAO"):
    out = [float(x) for x in res.get_BAO(Z)[:, 1]]
elif mode == "get_BAO_resParams":
    out = [float(x) for x in res.get_BAO(Z, res.Params)[:, 1]]
elif mode == "get_BAO_userparams":
    out = [float(x) for x in res.get_BAO(Z, cp)[:, 1]]
elif mode == "get_background_resParams":
    out = H(camb.get_background(res.Params))
elif mode == "get_background_resParams_copy":
    out = H(camb.get_background(res.Params.copy()))
elif mode == "calc_background_self":
    res.calc_background(res.Params)
    out = H(res)
elif mode == "user_cp_twice":
    out = H(camb.get_background(cp))
elif mode == "set_params_then_recalc":
    r2 = camb.CAMBdata()
    r2.set_params(res.Params)
    r2.calc_background(cp)
    out = H(r2)
print("RESULT", json.dumps(out), flush=True)
"""

MODES = [
    "fresh",
    "get_BAO_default",
    "get_BAO_resParams",
    "get_BAO_userparams",
    "get_background_resParams",
    "get_background_resParams_copy",
    "calc_background_self",
    "user_cp_twice",
    "set_params_then_recalc",
]


def tree_path(t: str) -> str:
    return t if os.path.isabs(t) else str(ROOT / "trees" / t)


def make_ini(model: str, tmp: str) -> str:
    ini = os.path.join(tmp, f"{model}.ini")
    Path(ini).write_text(
        f"DEFAULT({PARAMS / (model + '.ini')})\noutput_root = {tmp}/{model}\n"
        "EFTCAMB_write_background = F\nfeedback_level = 0\nEFTCAMB_feedback_level = 0\n"
    )
    return ini


def run(tree: str, ini: str, mode: str, timeout: int = 900) -> dict:
    env = dict(os.environ, OMP_NUM_THREADS="2")
    env.pop("PYTHONPATH", None)
    try:
        p = subprocess.run(
            [sys.executable, "-c", CHILD, tree, ini, mode, json.dumps(Z), str(Path(__file__).resolve().parents[1])],
            capture_output=True, text=True, env=env, cwd=os.path.dirname(ini), timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return {"mode": mode, "rc": "timeout"}
    lines = (p.stdout + p.stderr).splitlines()
    res = [ln for ln in lines if ln.startswith("RESULT")]
    flag = [ln.split()[1] for ln in lines if ln.startswith("EFTFLAG")]
    err = [ln.strip() for ln in lines if "Fortran runtime error" in ln or ln.startswith("At line")]
    pyerr = [ln for ln in lines if "Error" in ln and "Fortran" not in ln][-1:]
    return {
        "mode": mode, "rc": p.returncode, "eftflag": flag[0] if flag else None,
        "H": json.loads(res[0][7:]) if res else None, "err": " | ".join(err) or " ".join(pyerr),
    }


def scenarios(a) -> None:
    tree = tree_path(a.tree)
    tmp = tempfile.mkdtemp(prefix="vrealloc_", dir=OUT)
    for model in a.ini:
        ini = make_ini(model, tmp)
        modes = MODES + (["full_get_BAO"] if a.full else [])
        if Path(tree).name == "B1":
            modes = [m for m in modes if m != "get_BAO_default"]
        with ThreadPoolExecutor(a.jobs) as ex:
            rows = list(ex.map(lambda m: run(tree, ini, m), modes))
        print(f"== tree {Path(tree).name}  model {model}  EFTflag={rows[0]['eftflag']}")
        ref = rows[0]["H"]
        for r in rows:
            same = "" if r["H"] is None or ref is None else f" maxrel={max(abs(x / y - 1) for x, y in zip(r['H'], ref)):.1e}"
            print(f"  {r['mode']:32s} rc={r['rc']!s:4s} H={r['H']}{same} {r['err'][:160]}")
        (OUT / f"scen_{Path(tree).name}_{model}.json").write_text(json.dumps(rows, indent=1))


def allmodels(a) -> None:
    tree = tree_path(a.tree)
    tmp = tempfile.mkdtemp(prefix="vrealloc_all_", dir=OUT)
    models = sorted(p.stem for p in PARAMS.glob("*.ini") if p.stem not in ("base_params", "hdsk_base_params"))

    def one(m):
        ini = make_ini(m, tmp)
        return m, run(tree, ini, "fresh"), run(tree, ini, a.mode)

    with ThreadPoolExecutor(a.jobs) as ex:
        res = list(ex.map(one, models))
    n_crash = 0
    summary = {}
    for m, f, r in res:
        crash = r["rc"] != 0
        n_crash += crash
        summary[m] = {"fresh": f, a.mode: r}
        print(f"{m:34s} EFTflag={r['eftflag']} fresh_rc={f['rc']} {a.mode}_rc={r['rc']} {r['err'][:110]}")
    print(f"TOTAL models={len(res)} {a.mode} nonzero rc={n_crash}")
    (OUT / f"allmodels_{Path(tree).name}_{a.mode}.json").write_text(json.dumps(summary, indent=1))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("scenarios")
    s.add_argument("--tree", default="A0")
    s.add_argument("--ini", nargs="+", default=["4_Des_F_R_LCDM_1"])
    s.add_argument("--jobs", type=int, default=4)
    s.add_argument("--full", action="store_true")
    m = sub.add_parser("allmodels")
    m.add_argument("--tree", default="A0")
    m.add_argument("--mode", default="get_BAO_default")
    m.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    {"scenarios": scenarios, "allmodels": allmodels}[a.cmd](a)


if __name__ == "__main__":
    main()
