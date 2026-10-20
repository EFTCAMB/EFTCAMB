"""Verify: Horndeski V_scf (designer_V=F) reads never-written heap. Uses macOS MallocPreScribble (fills new
malloc blocks with 0xAA) as an independent way to expose uninitialised reads, without rebuilding."""
import sys
import numpy as np
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness

OUT = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/uninit"
cheap = {"get_transfer": "F", "get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "100"}
for model in ("5_hdsk_TG_1", "5_hdsk_TG_3"):
    ini = harness.PARAMS_DIR / f"{model}.ini"
    res = {}
    for tag, env in (("plain", {}), ("prescribble", {"MallocPreScribble": "1"}), ("scribble", {"MallocScribble": "1", "MallocPreScribble": "1"})):
        npz = f"{OUT}/{model}_{tag}.npz"
        st = harness.run_model("A0", ini, npz, overrides=cheap, threads=1, extra_env=env)
        d = harness.load(npz)
        v = np.asarray(d["eft_V_scf"])
        res[tag] = v
        print(model, tag, st["status"], "V_scf: nonzero", int(np.count_nonzero(v)), "/", v.size,
              "unique", np.unique(v)[:4], "phi_scf[-1]", d["eft_phi_scf"][-1], "TT[50]", d["cl_total"][50,0])
