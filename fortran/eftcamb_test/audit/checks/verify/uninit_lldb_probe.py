"""lldb script: at the first EFTCAMBTimestepCacheIsNan call of every EFTCAMB_Stability_Check, print the raw contents of
eft_cache%EFTcdotdot (offset 240) and eft_cache%EFTGamma1PP (offset 288) as seen by the NaN check in a release build.

Loaded by verify_uninit_lldb.py with `command script import`. The cache is a class(TEFTCAMB_timestep_cache) dummy, so
x0 points to the class container whose first word is the data address; field `a` (offset 0) is printed as a sanity check.
"""

import struct

import lldb

STATE = {"armed": False, "check_bp": None, "nan_bp": None}


def __lldb_init_module(debugger, _dict):
    target = debugger.GetSelectedTarget()
    b1 = target.BreakpointCreateByName("__eftcamb_stability_MOD_eftcamb_stability_check")
    b1.SetScriptCallbackFunction("uninit_lldb_probe.on_check")
    b2 = target.BreakpointCreateByName("__eftcamb_cache_MOD_eftcambtimestepcacheisnan")
    b2.SetScriptCallbackFunction("uninit_lldb_probe.on_isnan")
    STATE["nan_bp"] = b2.GetID()


def on_check(frame, bp_loc, _dict):
    STATE["armed"] = True
    return False


def on_isnan(frame, bp_loc, _dict):
    if not STATE["armed"]:
        return False
    STATE["armed"] = False
    proc = frame.GetThread().GetProcess()
    err = lldb.SBError()
    x0 = frame.FindRegister("x0").GetValueAsUnsigned()
    data = struct.unpack("<Q", proc.ReadMemory(x0, 8, err))[0]
    vals = {}
    for off, name in ((0, "a"), (240, "EFTcdotdot"), (288, "EFTGamma1PP")):
        raw = proc.ReadMemory(data + off, 8, err)
        vals[name] = (struct.unpack("<d", raw)[0], raw[::-1].hex())
    print("PROBE " + " ".join(f"{k}={v[0]!r}(0x{v[1]})" for k, v in vals.items()), flush=True)
    return False
