"""t3-gr: source-window (number counts / galaxy lensing) Cls of EFT models at their GR point vs LCDM.

Base params have num_redshiftwindows=2 (counts at z=0.5 with redshift distortions, lensing window).
The harness does not save source-window Cls, so this script has its own worker.

python a7_windows.py [case,...]
python a7_windows.py _worker TREE INI OUT.npz JSON_PYTHON_SET
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def worker(tree: str, ini: str, out: str, pset: str) -> None:
    sys.path.insert(0, tree)
    import numpy as np

    import camb

    assert Path(camb.__file__).resolve().is_relative_to(Path(tree).resolve())
    cp = camb.read_ini(ini, no_validate=True)
    for k, v in json.loads(pset).items():
        obj = cp
        parts = k.split(".")
        for p_ in parts[:-1]:
            obj = getattr(obj, p_)
        setattr(obj, parts[-1], v)
    r = camb.get_results(cp)
    d = r.get_source_cls_dict(raw_cl=True, lmax=1000)
    np.savez(out, **{k: np.asarray(v) for k, v in d.items()},
             EFTflag=np.array(getattr(getattr(cp, "EFTCAMB", None), "EFTflag", -1)))


def main() -> None:
    import numpy as np

    from common import GR_INI, PARAMS, WORK, harness

    sub = WORK / "runs" / "a7"
    sub.mkdir(parents=True, exist_ok=True)
    pe = PARAMS / "2_PEFT_Omega_const_1.ini"

    def wrap(name, base_ini, ov):
        p = sub / f"{name}.ini"
        p.write_text("\n".join([f"DEFAULT({base_ini})", f"output_root = {sub / name}",
                                "EFTCAMB_write_background = F", "feedback_level = 0"]
                               + [f"{k} = {v}" for k, v in ov.items()]) + "\n")
        return p

    cases = {
        "B0_GR": ("B0", GR_INI, {}, {}),
        "A0_GR": ("A0", GR_INI, {}, {"Reion.max_redshift": 50.0}),
        "A0_PEFTzero": ("A0", pe, {"PureEFTmodelOmega": "0"}, {"Reion.max_redshift": 50.0}),
        "A0_PEFTzero_forced": ("A0", pe, {"PureEFTmodelOmega": "0", "EFTCAMB_GR_threshold": "-1"},
                               {"Reion.max_redshift": 50.0}),
        "A0_PEFT_Om1e-7": ("A0", pe, {"EFTOmega0": "1e-7"}, {"Reion.max_redshift": 50.0}),
        "A0_PEFT_Om1e-9": ("A0", pe, {"EFTOmega0": "1e-9"}, {"Reion.max_redshift": 50.0}),
        "A0_PEFTzero_noRSD": ("A0", pe, {"PureEFTmodelOmega": "0", "counts_redshift": "F"},
                              {"Reion.max_redshift": 50.0}),
        "B0_GR_noRSD": ("B0", GR_INI, {"counts_redshift": "F"}, {}),
        "A0_fR_B0_1e-9": ("A0", PARAMS / "4_Des_F_R_LCDM_1.ini", {"EFTB0": "1e-9"}, {"Reion.max_redshift": 50.0}),
        "A0_fR_B0_1e-6": ("A0", PARAMS / "4_Des_F_R_LCDM_1.ini", {"EFTB0": "1e-6"}, {"Reion.max_redshift": 50.0}),
    }
    if len(sys.argv) > 1:
        cases = {k: cases[k] for k in sys.argv[1].split(",")}

    def one(item):
        name, (tree, ini, ov, ps) = item
        out = sub / f"{name}.npz"
        if out.exists():
            return name, 0, ""
        env = dict(os.environ, OMP_NUM_THREADS="2")
        env.pop("PYTHONPATH", None)
        p = subprocess.run([sys.executable, __file__, "_worker", str(harness.tree_path(tree)), str(wrap(name, ini, ov)),
                            str(out), json.dumps(ps)], env=env, capture_output=True, text=True, cwd=sub)
        return name, p.returncode, (p.stdout + p.stderr)[-1500:]

    with ThreadPoolExecutor(int(os.environ.get("T3_JOBS", "4"))) as ex:
        for name, rc, log in ex.map(one, cases.items()):
            if rc:
                print(name, "FAILED", log)

    def load(n):
        with np.load(sub / f"{n}.npz") as f:
            return {k: f[k] for k in f.files}

    ref = load("B0_GR")
    keys = [k for k in ref if k not in ("EFTflag", "PxP")]
    for n in cases:
        if not (sub / f"{n}.npz").exists():
            continue
        refn = load("B0_GR_noRSD") if n.endswith("noRSD") else ref
        d = load(n)
        devs = {}
        for k in keys:
            a, b = d[k][2:], refn[k][2:]
            norm = np.max(np.abs(b))
            devs[k] = float(np.max(np.abs(a - b)) / norm) if norm else 0.0
        print(f"{n:22s} EFTflag={int(d['EFTflag']):2d} " + " ".join(f"{k}={v:.1e}" for k, v in devs.items()))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "_worker":
        worker(*sys.argv[2:6])
    else:
        main()
