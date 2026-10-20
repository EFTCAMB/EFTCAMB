"""Independent check of ic23-stale-cache-gammas on the unmodified A0 tree: with EFT_IC_type=3 the ICs read
EV%eft_cache Gammas left by whatever mode the thread ran before, so results should depend on thread count."""
import json, sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import harness
A = harness.ROOT if hasattr(harness, "ROOT") else None
S = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
tree = S + "/trees/A0"
out = pathlib.Path(S + "/work/verify/ic23"); out.mkdir(parents=True, exist_ok=True)
P = harness.REPO / "fortran/eftcamb_test/parameters"
models = sys.argv[1:]
for m in models:
    res = {}
    for tag, thr, ic in (("ic3_t1", 1, 3), ("ic3_t1b", 1, 3), ("ic3_t2", 2, 3), ("ic1_t1", 1, 1), ("ic1_t2", 2, 1)):
        f = out / f"{m}_on1e-8_{tag}.npz"
        st = harness.run_model(tree, str(P / f"{m}.ini"), str(f), overrides={"EFT_IC_type": ic, "EFTCAMB_turn_on_time": "1e-8"}, threads=thr)
        res[tag] = st["status"]
    line = {"model": m, "status": res}
    def cmp(a, b):
        if res[a] != "ok" or res[b] != "ok":
            return None
        c = harness.compare(harness.load(str(out / f"{m}_on1e-8_{a}.npz")), harness.load(str(out / f"{m}_on1e-8_{b}.npz")))
        return {"worst_cl": c["worst_cl"], "worst_pk": c["worst_pk"], "bit": c["bit_identical"]}
    line["ic3 t1 vs t1b"] = cmp("ic3_t1", "ic3_t1b")
    line["ic3 t1 vs t2"] = cmp("ic3_t1", "ic3_t2")
    line["ic1 t1 vs t2"] = cmp("ic1_t1", "ic1_t2")
    line["ic3 vs ic1 (t1)"] = cmp("ic3_t1", "ic1_t1")
    print(json.dumps(line, default=str), flush=True)
