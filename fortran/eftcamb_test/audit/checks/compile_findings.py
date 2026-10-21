"""Compile fortran/eftcamb_test/audit/FINDINGS.md from the audit workflow journals.

Reads every {"type": "result"} record: track/physics results (verdict, checks, findings, matrix) and verifier verdicts
(labels ``verify:<owner>:<slug>:<lens>``), merges them with the orchestrator's own findings, and writes the ledger.

Usage: python compile_findings.py JOURNAL.jsonl [JOURNAL2.jsonl ...] [--out FINDINGS.md]
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]

# findings made directly by the orchestrator while setting up the audit
ORCHESTRATOR = [
    {
        "slug": "python-string-getters-hidden-length",
        "title": "Python calls character(len=200) getters without gfortran's hidden length argument",
        "severity": "S3",
        "file": "fortran/eftcamb/11_EFTCAMB_python.f90",
        "line": 248,
        "observed": "-fcheck build: 'Actual string length is shorter than the declared one for dummy argument "
        "model_name (127/200)' on every EFTCAMB.model_name() call; same for get_read_parameters, param names/labels, "
        "set_model_name.",
        "expected": "Explicit-length interface (bind(c) with a length argument, or the length passed from ctypes).",
        "root_cause": "ctypes passes only the buffer pointer; gfortran reads the hidden length from a garbage register. "
        "Release builds assume len=200 and work only because Python allocates exactly 200 bytes.",
        "proposed_fix": "bind(c) wrappers taking (buffer, length) and copying with explicit length; or pass "
        "c_long(200) as the trailing hidden argument from camb/eftcamb.py.",
        "fix_class": "safe-local",
        "confidence": "high",
        "track": "orchestrator",
        "verdicts": [],
    },
    {
        "slug": "upstream-thermo-size-unallocated",
        "title": "upstream: size() of an unallocated array in ThermoData init (non-short-circuit .and.)",
        "severity": "S3",
        "file": "fortran/results.f90",
        "line": 2168,
        "observed": "-fcheck build aborts every run: 'Allocatable argument this is not allocated' at "
        "`if (allocated(this%tb) .and. this%nthermo /= size(this%tb))`.",
        "expected": "Nested if; Fortran does not guarantee short-circuit evaluation of .and.",
        "root_cause": "Upstream CAMB code (same line in branch camb, line 2035).",
        "proposed_fix": "if (allocated(this%tb)) then; if (this%nthermo /= size(this%tb)) then ... (report upstream).",
        "fix_class": "safe-local",
        "confidence": "high",
        "track": "orchestrator",
        "verdicts": [],
    },
    {
        "slug": "eftflag0-read-parameters-only-eftflag",
        "title": "With EFTflag=0, read_parameters() returns only EFTflag, so set_params(EFTflag=0, <other EFT key>) raises "
        "CAMBUnknownArgumentError",
        "severity": "S3",
        "file": "camb/camb.py",
        "line": 224,
        "observed": "set_params(EFTflag=0, EFT_additional_priors=...) raises CAMBUnknownArgumentError (reported by the "
        "user, 2026-10-19); the pre-rebase code hid it through a stale global parameter cache.",
        "expected": "EFT keys accepted (and ignored) when EFTflag=0, as they are by the ini path.",
        "root_cause": "set_params marks as used only the keys returned by read_parameters, which lists only EFTflag "
        "when no model is allocated.",
        "proposed_fix": "In EFT_do_set, also mark as used every key the EFTCAMB ini reader recognises (or every key "
        "starting with EFT/known model prefixes) when EFTflag=0.",
        "fix_class": "safe-local",
        "confidence": "medium",
        "track": "user-report",
        "verdicts": [],
    },
]

SEV_ORDER = {"S1": 0, "S2": 1, "S3": 2}


def load_records(journals: list[Path]) -> tuple[dict, dict]:
    """Return {label: result} for the last result of each label, and {label: agentId}."""
    results: dict[str, object] = {}
    for j in journals:
        labels: dict[str, str] = {}
        for line in j.read_text().splitlines():
            rec = json.loads(line)
            if rec["type"] == "started":
                labels[rec["key"]] = rec["label"]
            elif rec["type"] == "result" and rec.get("result") is not None:
                results[labels.get(rec["key"], rec["key"])] = rec["result"]
    return results, {}


def status_of(f: dict) -> str:
    vs = f.get("verdicts") or []
    if f.get("track") == "orchestrator":
        return "confirmed by orchestrator (reproduced)"
    if not vs:
        return "unverified" if f["severity"] == "S3" else "verification pending"
    kinds = [v["verdict"] for v in vs]
    if "REFUTED" in kinds and "CONFIRMED" not in kinds:
        return "refuted"
    if all(k == "CONFIRMED" for k in kinds):
        return "confirmed"
    return "/".join(kinds).lower()


def effective_severity(f: dict) -> str:
    sev = [v.get("corrected_severity") for v in f.get("verdicts") or [] if v.get("corrected_severity")]
    if not sev:
        return f["severity"]
    if all(s == "not-a-bug" for s in sev):
        return "not-a-bug"
    known = [s for s in sev if s in SEV_ORDER]
    return min(known, key=SEV_ORDER.get) if known else f["severity"]


def md_escape(s: object) -> str:
    return str(s).replace("|", "\\|").replace("\n", " ")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("journals", nargs="+", type=Path)
    p.add_argument("--out", type=Path, default=REPO / "fortran" / "eftcamb_test" / "audit" / "FINDINGS.md")
    a = p.parse_args()
    results, _ = load_records(a.journals)

    tracks = {k: v for k, v in results.items() if isinstance(v, dict) and ("track" in v or "matrix" in v)}
    verdicts: dict[tuple[str, str], list] = defaultdict(list)
    for label, v in results.items():
        if label.startswith("verify:") and isinstance(v, dict):
            _, owner, slug, lens = label.split(":", 3)
            verdicts[(owner, slug)].append({**v, "lens": lens})

    findings = list(ORCHESTRATOR)
    for label, r in tracks.items():
        owner = r.get("track") or r.get("id") or label
        for f in r.get("findings") or []:
            findings.append({**f, "track": owner, "verdicts": verdicts.get((owner, f["slug"]), [])})
    for i, f in enumerate(sorted(findings, key=lambda f: (SEV_ORDER.get(effective_severity(f), 3), f["track"]))):
        f["id"] = f"F-{i + 1:03d}"
    findings.sort(key=lambda f: f["id"])

    out: list[str] = ["# EFTCAMB audit findings", ""]
    out.append(
        "Code under audit: branch `eftcamb-audit` at 5bef49a (EFTCAMB on CAMB v2). Generated by "
        "`fortran/eftcamb_test/audit/checks/compile_findings.py` from the audit workflow results; every S1/S2 "
        "finding has a reproducer under `fortran/eftcamb_test/audit/repro/`."
    )
    out.append("")
    counts: dict[str, int] = defaultdict(int)
    for f in findings:
        if status_of(f) != "refuted":
            counts[effective_severity(f)] += 1
    out.append(
        "**Counts (excluding refuted):** "
        + ", ".join(f"{s}: {counts.get(s, 0)}" for s in ("S1", "S2", "S3", "not-a-bug"))
    )
    out.append("")
    out.append("## Track verdicts")
    out.append("")
    out.append("| Track | Verdict | Findings | Not covered |")
    out.append("|---|---|---|---|")
    for label, r in sorted(tracks.items()):
        owner = r.get("track") or r.get("id") or label
        nc = "; ".join(r.get("not_covered") or []) or "-"
        out.append(f"| {owner} | {r.get('verdict')} | {len(r.get('findings') or [])} | {md_escape(nc)[:400]} |")
    out.append("")
    out.append("## Ledger")
    out.append("")
    out.append("| ID | Sev | Status | Title | Location | Track | Fix class |")
    out.append("|---|---|---|---|---|---|---|")
    for f in findings:
        loc = f"{f['file']}:{f.get('line', '')}".rstrip(":")
        out.append(
            f"| {f['id']} | {effective_severity(f)} | {status_of(f)} | {md_escape(f['title'])} | `{loc}` | "
            f"{f['track']} | {f.get('fix_class', '')} |"
        )
    out.append("")
    out.append("## Details")
    for f in findings:
        out.append("")
        out.append(f"### {f['id']} [{effective_severity(f)}, {status_of(f)}] {f['title']}")
        out.append("")
        loc = f"{f['file']}:{f.get('line', '')}".rstrip(":")
        out.append(f"- **Location:** `{loc}`" + (f" (also: {', '.join(f['other_locations'])})" if f.get("other_locations") else ""))
        out.append(f"- **Track:** {f['track']}; reported severity {f['severity']}, confidence {f.get('confidence', '')}")
        if f.get("repro_script") or f.get("repro_command"):
            out.append(f"- **Reproduce:** `{f.get('repro_command') or f.get('repro_script')}`")
        out.append(f"- **Observed:** {f['observed']}")
        out.append(f"- **Expected:** {f['expected']}")
        if f.get("code_expression"):
            out.append(f"- **Code:** {f['code_expression']}")
        if f.get("derived_expression"):
            out.append(f"- **Derived:** {f['derived_expression']}")
        out.append(f"- **Root cause:** {f['root_cause']}")
        out.append(f"- **Proposed fix ({f.get('fix_class', '')}):** {f['proposed_fix']}")
        for v in f.get("verdicts") or []:
            out.append(
                f"- **Verifier ({v['lens']}):** {v['verdict']}, severity {v.get('corrected_severity')}; {v['reason']}"
            )
    out.append("")
    out.append("## Checks per track")
    for label, r in sorted(tracks.items()):
        owner = r.get("track") or r.get("id") or label
        out.append("")
        out.append(f"### {owner}: {r.get('verdict')}")
        out.append("")
        out.append(md_escape(r.get("summary", "")))
        out.append("")
        if r.get("checks"):
            out.append("| Check | Value | Threshold | Pass |")
            out.append("|---|---|---|---|")
            for c in r["checks"]:
                out.append(
                    f"| {md_escape(c['name'])} | {md_escape(c['value'])} | {md_escape(c.get('threshold', ''))} | "
                    f"{'PASS' if c['pass'] else 'FAIL'} |"
                )
        if r.get("commands"):
            out.append("")
            out.append("Commands:")
            out.extend(f"- `{md_escape(c)}`" for c in r["commands"])
    matrix = [dict(m, owner=r.get("id")) for r in tracks.values() for m in r.get("matrix") or []]
    if matrix:
        out.append("")
        out.append("## Physics matrix")
        out.append("")
        out.append("| Model | Check | Value | Threshold | Pass | Source |")
        out.append("|---|---|---|---|---|---|")
        for m in sorted(matrix, key=lambda m: (m["model"], m["check"])):
            out.append(
                f"| {m['model']} | {md_escape(m['check'])} | {md_escape(m['value'])} | {md_escape(m.get('threshold', ''))} "
                f"| {'PASS' if m['pass'] else 'FAIL'} | {m['owner']} |"
            )
    a.out.write_text("\n".join(out) + "\n")
    print(f"wrote {a.out}: {len(findings)} findings, {len(tracks)} track results, {len(matrix)} matrix entries")
    print(dict(counts))


if __name__ == "__main__":
    main()
