"""Condense interrupted subagent transcripts into progress notes an agent can reread when resuming.

Keeps the agent's own text, the commands/files it ran or wrote (truncated), and short tool-result tails.

Usage: python condense_transcript.py MAP.json OUTDIR   (MAP: {label: transcript.jsonl})
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


def condense(path: Path, max_result: int = 600) -> str:
    out: list[str] = []
    for line in path.read_text(errors="replace").splitlines():
        try:
            rec = json.loads(line)
        except json.JSONDecodeError:
            continue
        msg = rec.get("message") or {}
        content = msg.get("content")
        if isinstance(content, str):
            content = [{"type": "text", "text": content}]
        for block in content or []:
            kind = block.get("type")
            if kind == "text" and msg.get("role") == "assistant":
                out.append(f"\n[note] {block['text'].strip()}")
            elif kind == "tool_use":
                inp = block.get("input", {})
                desc = inp.get("command") or inp.get("file_path") or json.dumps(inp)[:300]
                out.append(f"[{block.get('name')}] {str(desc)[:500]}")
            elif kind == "tool_result":
                res = block.get("content")
                if isinstance(res, list):
                    res = " ".join(c.get("text", "") for c in res if isinstance(c, dict))
                res = str(res or "").strip()
                if res:
                    out.append(f"  -> {res[-max_result:]}")
    return "\n".join(out)


def main() -> None:
    mapping = json.loads(Path(sys.argv[1]).read_text())
    outdir = Path(sys.argv[2])
    outdir.mkdir(parents=True, exist_ok=True)
    for label, transcript in mapping.items():
        p = Path(transcript)
        if not p.exists() or label.startswith("verify"):
            continue
        text = condense(p)
        dest = outdir / (label.replace(":", "_").replace("/", "_") + ".md")
        dest.write_text(text)
        print(f"{label:30s} {len(text) // 1024:6d} KB -> {dest}")


if __name__ == "__main__":
    main()
