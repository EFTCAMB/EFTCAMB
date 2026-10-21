"""t3-gr: is the curvature problem of EFT models inherited from the pre-rebase EFTCAMB (B1)?"""
from common import PARAMS, harness, run, fmt

ov = {"PureEFTmodelOmega": "0", "omk": "0.05", "num_redshiftwindows": "0"}
p1 = run("B1", "B1_PEFTzero_omk", ini=PARAMS / "2_PEFT_Omega_const_1.ini", overrides=ov, sub="b3")
p0 = run("B1", "B1_GR_omk", ini=PARAMS / "1_EFT_GR.ini", overrides={"omk": "0.05", "num_redshiftwindows": "0"},
         sub="b3")
r = harness.compare(p1, p0)
print("B1 PEFT(zero) omk=0.05 vs B1 GR omk=0.05:", fmt(r))
a, b = harness.load(p1), harness.load(p0)
for l in (2, 10, 30, 100):
    print(l, a["cl_unlensed_scalar"][l, 0] / b["cl_unlensed_scalar"][l, 0] - 1)
