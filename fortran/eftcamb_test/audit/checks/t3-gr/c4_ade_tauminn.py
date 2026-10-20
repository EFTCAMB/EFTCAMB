"""t3-gr: 5_ADE_1 (EFTCAMB_turn_on_time = 1e-10) fails with 'EFTCAMB starts before thermo tauminn' when CMB settings
are reduced. Which setting moves ThermoData%tauminn past tau(a=1e-10)? One override at a time; also EFTflag=0 control."""
import re
from concurrent.futures import ThreadPoolExecutor

from common import PARAMS, RUNS, harness

INI = PARAMS / "5_ADE_1.ini"
CASES = {"default": {}, "nolens": {"do_lensing": "F"}, "notensor": {"get_tensor_cls": "F"},
         "lmax2000": {"l_max_scalar": "2000"}, "lmax2500_lens": {"l_max_scalar": "2500"},
         "notransfer": {"get_transfer": "F"}}


def one(item):
    n, ov = item
    out = RUNS / "c4" / f"ADE_{n}.npz"
    if not out.with_suffix(".json").exists():
        harness.run_model("A0", INI, out, overrides={**ov, "get_tensor_cls": ov.get("get_tensor_cls", "T")},
                          threads=1, timeout=900)
    st = harness.json.loads(out.with_suffix(".json").read_text())
    log = out.with_suffix(".log").read_text(errors="replace")
    m = re.findall(r"tauminn = ([^\n]*)", log)
    return n, st["status"], m[:1]


with ThreadPoolExecutor(3) as ex:
    for r in ex.map(one, CASES.items()):
        print(r)
