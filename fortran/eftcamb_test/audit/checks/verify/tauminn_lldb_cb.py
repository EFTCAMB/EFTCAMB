"""lldb callbacks for tauminn_lldb.py: print double arguments of Thermo_Init (taumin) and GetTauStart (q)."""

import lldb


def _dbl(frame, reg):
    proc = frame.GetThread().GetProcess()
    addr = frame.FindRegister(reg).GetValueAsUnsigned()
    err = lldb.SBError()
    data = proc.ReadMemory(addr, 8, err)
    import struct

    return struct.unpack("<d", data)[0] if err.Success() else float("nan")


def thermo_cb(frame, bp_loc, extra_args, internal_dict):
    caller = frame.GetThread().GetFrameAtIndex(1).GetFunctionName()
    print(f"AUDIT Thermo_Init taumin={_dbl(frame, 'x2'):.6e} tauminn=0.95*taumin={0.95 * _dbl(frame, 'x2'):.6e}"
          f" caller={caller}", flush=True)
    return False


def taustart_cb(frame, bp_loc, extra_args, internal_dict):
    caller = frame.GetThread().GetFrameAtIndex(1).GetFunctionName()
    print(f"AUDIT GetTauStart q={_dbl(frame, 'x0'):.6e} caller={caller}", flush=True)
    return False
