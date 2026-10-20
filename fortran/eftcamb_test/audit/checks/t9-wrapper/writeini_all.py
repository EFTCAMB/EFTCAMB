"""Track 9 (iv): CAMBparams.write_ini with EFTCAMB on, for every test ini (no Boltzmann runs).

Per model: p = read_ini(model); d1 = p.EFTCAMB.read_parameters(); write_ini(tmp, validate=True) (records the
exception); p2 = read_ini(tmp) (written with validate=False); d2 = p2.EFTCAMB.read_parameters().
Reports: number of EFT keys in the written file, validate outcome, EFTflag after re-read, d2 == d1 (output_root
excluded).  Usage: python writeini_all.py [tree]   (default A0)
"""

import json
import sys

from t9lib import run_snippet

CODE = r"""
import tempfile, glob
res = {}
for ini in sorted(glob.glob(str(PARAMS / '*.ini'))):
    stem = Path(ini).stem
    if 'base_params' in stem:
        continue
    p = camb.read_ini(ini)
    d1 = {k: v for k, v in p.EFTCAMB.read_parameters().items() if k != 'output_root'}
    with tempfile.TemporaryDirectory() as tmp:
        f = os.path.join(tmp, 'w.ini')
        try:
            p.write_ini(f, validate=True)
            val = 'ok'
        except Exception as e:
            val = type(e).__name__ + ': ' + str(e).split('):')[-1][:80]
        p.write_ini(f, validate=False)
        keys = [ln.split('=')[0].strip() for ln in open(f) if '=' in ln]
        n_eft = sum(k in d1 for k in keys)
        p2 = camb.read_ini(f)
        d2 = {k: v for k, v in p2.EFTCAMB.read_parameters().items() if k != 'output_root'}
    res[stem] = dict(EFTflag=int(p.EFTCAMB.EFTflag), n_eft_read=len(d1), n_eft_written=n_eft, validate=val,
                     EFTflag_reread=int(p2.EFTCAMB.EFTflag), same=(d1 == d2))
emit(res=res)
"""

if __name__ == "__main__":
    tree = sys.argv[1] if len(sys.argv) > 1 else "A0"
    r = run_snippet(CODE, tree_name=tree, threads=1, timeout=1500)
    res = r["result"].get("res", {})
    if r["status"] != "ok":
        print(r["status"], "\n".join(r["tail"][-8:]))
    eft = {k: v for k, v in res.items() if v["EFTflag"] != 0}
    print(f"{tree}: models={len(res)} (EFT on: {len(eft)}); round trip identical: "
          f"{sum(v['same'] for v in eft.values())}/{len(eft)}; validate ok: {sum(v['validate'] == 'ok' for v in eft.values())}"
          f"/{len(eft)}; EFT keys written (total) {sum(v['n_eft_written'] for v in eft.values())} of "
          f"{sum(v['n_eft_read'] for v in eft.values())}")
    for k, v in res.items():
        if not v["same"] or v["validate"] != "ok":
            print(f"  {k:30s} {json.dumps(v)[:200]}")
