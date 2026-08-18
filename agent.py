#!/usr/bin/env python3
"""
ClawBio Conversational Agent
=============================
A terminal agent powered by any OpenAI-compatible LLM (Nebius, OpenAI,
Ollama, Groq, etc.) that can call ClawBio skills as tools.

The agent reads a natural language prompt, decides which skills to run,
executes them, and synthesises a biological answer from the results.

Usage:
    python agent.py "Analyse my human genome and find if these genes share structural homology"
    python agent.py --interactive          # chat loop
    python agent.py --demo                 # run with a canned example prompt

Config (set in .env or environment):
    LLM_API_KEY     your Nebius / OpenAI / Groq key
    LLM_BASE_URL    https://api.studio.nebius.com/v1  (or any OpenAI-compat base)
    CLAWBIO_MODEL   meta-llama/Meta-Llama-3.1-70B-Instruct  (or any model name)
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import textwrap
from datetime import datetime
from pathlib import Path

# ---------------------------------------------------------------------------
# Bootstrap
# ---------------------------------------------------------------------------

_ROOT = Path(__file__).resolve().parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

try:
    from dotenv import load_dotenv
    load_dotenv(_ROOT / ".env")
    load_dotenv()
except ImportError:
    pass  # dotenv optional — env vars can be set directly

try:
    from openai import OpenAI
except ImportError:
    print("openai package not found. Run:  pip install openai")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

LLM_API_KEY   = (
    os.environ.get("LLM_API_KEY")
    or os.environ.get("OPENAI_API_KEY")
    or os.environ.get("NEBIUS_API_KEY")
    or ""
)
LLM_BASE_URL  = os.environ.get("LLM_BASE_URL", "")
MODEL         = os.environ.get("CLAWBIO_MODEL", "meta-llama/Llama-3.3-70B-Instruct")

# Auto-detect Nebius base URL if only NEBIUS_API_KEY is set
if not LLM_BASE_URL and os.environ.get("NEBIUS_API_KEY") and not os.environ.get("LLM_API_KEY"):
    LLM_BASE_URL = "https://api.studio.nebius.com/v1"
SKILLS_DIR    = _ROOT / "skills"
OUTPUT_DIR    = _ROOT / "output" / f"agent_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

if not LLM_API_KEY:
    print(
        "LLM_API_KEY not set.\n"
        "For Nebius: export LLM_API_KEY=<your-key> LLM_BASE_URL=https://api.studio.nebius.com/v1\n"
        "For OpenAI: export LLM_API_KEY=<your-key>"
    )
    sys.exit(1)

client_kwargs: dict = {"api_key": LLM_API_KEY}
if LLM_BASE_URL:
    client_kwargs["base_url"] = LLM_BASE_URL

llm = OpenAI(**client_kwargs)

# ---------------------------------------------------------------------------
# System prompt
# ---------------------------------------------------------------------------

SYSTEM_PROMPT = textwrap.dedent("""
You are a ClawBio bioinformatics agent — an expert computational biologist
who reasons carefully about genomic and structural biology questions.

You have access to a set of ClawBio tools. When the user asks a biological
question, you:
1. Think about which tool(s) to call.
2. Call the tools (they run locally — no data leaves the machine).
3. Read the results carefully.
4. Synthesise a clear, evidence-based biological answer.

Rules:
- Always cite tool output when making claims. Do not invent biology.
- If a tool is unavailable (binary not installed), say so clearly and explain
  how to install it, then continue with what you can.
- Keep answers scientific but accessible.
- Include the ClawBio disclaimer at the end of any health-related answer:
  "ClawBio is a research and educational tool. It is not a medical device
   and does not provide clinical diagnoses. Consult a healthcare professional
   before making any medical decisions."
""").strip()

# ---------------------------------------------------------------------------
# Tool definitions (OpenAI function-calling format)
# ---------------------------------------------------------------------------

TOOLS = [
    {
        "type": "function",
        "function": {
            "name": "run_foldseek",
            "description": (
                "Search a protein structure (CIF or PDB) against major structural "
                "databases (PDB, AlphaFold DB) and rank hits by TM-score. "
                "Use this to answer questions like: 'do these proteins share structural "
                "homology?', 'is this protein fold known?', 'find structurally similar proteins'. "
                "Run this BEFORE struct_predict when you want to check if a structure already "
                "exists — if TM-score >= 0.9 prediction is unnecessary."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "input_cif": {
                        "type": "string",
                        "description": "Path to a CIF or PDB file. Leave empty to use demo data.",
                    },
                    "databases": {
                        "type": "string",
                        "description": "Comma-separated databases: pdb, afdb, esm. Default: pdb",
                    },
                    "min_tmscore": {
                        "type": "number",
                        "description": "Minimum TM-score threshold for reported hits. Default: 0.3",
                    },
                },
                "required": [],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "struct_predict",
            "description": (
                "Predict the 3-D structure of a protein from a sequence using Boltz-2. "
                "Use this when the user provides a sequence and wants to know the structure, "
                "OR after run_foldseek confirms the structure is novel (TM-score < 0.9). "
                "Do NOT predict if foldseek already found a near-identical structure."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "input_yaml": {
                        "type": "string",
                        "description": "Path to Boltz-2 YAML input. Leave empty for demo.",
                    },
                },
                "required": [],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "run_pharmgx",
            "description": (
                "Run pharmacogenomics analysis on a genetic data file (23andMe, AncestryDNA, VCF). "
                "Reports gene-drug interactions for CYP2D6, CYP2C19, VKORC1, and more."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "input_file": {
                        "type": "string",
                        "description": "Path to genetic data file. Leave empty for demo.",
                    },
                },
                "required": [],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "run_gwas_lookup",
            "description": (
                "Look up a genetic variant (rsID) across 9 GWAS and annotation databases. "
                "Use for questions like: 'what does rs1234 do?', 'is this variant associated with disease?'"
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "rsid": {
                        "type": "string",
                        "description": "The rsID to look up (e.g. rs3798220).",
                    },
                },
                "required": ["rsid"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "run_skill",
            "description": (
                "Run any other registered ClawBio skill by name. "
                "Available: pharmgx, equity, nutrigx, prs, gwas, compare, scrna, rnaseq, "
                "phylo, analyze-fasta, foldseek, struct-predictor, metagenomics, clinpgx. "
                "Use this for skills not covered by the specific tools above."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "skill": {
                        "type": "string",
                        "description": "Skill name as registered in clawbio.py (e.g. 'phylo', 'prs')",
                    },
                    "input_file": {
                        "type": "string",
                        "description": "Path to input file. Leave empty for demo.",
                    },
                    "extra_args": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "Extra CLI arguments e.g. ['--trait', 'type 2 diabetes']",
                    },
                },
                "required": ["skill"],
            },
        },
    },
]

# ---------------------------------------------------------------------------
# Tool execution
# ---------------------------------------------------------------------------

CLAWBIO_PY = _ROOT / "clawbio.py"


def _run(cmd: list[str], label: str) -> dict:
    """Run a subprocess and return a structured result dict."""
    print(f"\n  ▶ {label}")
    print(f"    $ {' '.join(str(c) for c in cmd)}")
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(_ROOT),
            capture_output=True,
            text=True,
            timeout=600,
        )
    except FileNotFoundError as exc:
        return {
            "success": False,
            "error": str(exc),
            "stdout": "",
            "stderr": str(exc),
        }
    except subprocess.TimeoutExpired:
        return {
            "success": False,
            "error": "Timed out after 600s",
            "stdout": "",
            "stderr": "timeout",
        }

    ok = proc.returncode == 0
    if ok:
        print(f"    ✓ Done")
    else:
        print(f"    ✗ Exit {proc.returncode}")
        if proc.stderr:
            print(f"    stderr: {proc.stderr[:300]}")

    return {
        "success": ok,
        "exit_code": proc.returncode,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
    }


def _skill_output_dir(skill: str) -> Path:
    d = OUTPUT_DIR / skill
    d.mkdir(parents=True, exist_ok=True)
    return d


def _load_result_json(out_dir: Path) -> dict | None:
    rj = out_dir / "result.json"
    if rj.exists():
        try:
            return json.loads(rj.read_text())
        except Exception:
            pass
    return None


def _load_report_md(out_dir: Path) -> str:
    rmd = out_dir / "report.md"
    if rmd.exists():
        text = rmd.read_text()
        # Trim to first 3000 chars so the LLM context stays manageable
        return text[:3000] + ("\n...[truncated]" if len(text) > 3000 else "")
    return ""


def execute_tool(name: str, args: dict) -> str:
    """Dispatch a tool call and return a JSON string for the LLM."""

    # ── run_foldseek ───────────────────────────────────────────────────────
    if name == "run_foldseek":
        out_dir = _skill_output_dir("foldseek")
        cmd = [sys.executable, str(SKILLS_DIR / "struct-predictor-foldseek" / "struct_predictor_foldseek.py")]
        input_cif = args.get("input_cif", "")
        if input_cif:
            cmd += ["--input", input_cif]
        else:
            cmd += ["--demo"]
        if args.get("databases"):
            cmd += ["--databases", args["databases"]]
        if args.get("min_tmscore") is not None:
            cmd += ["--min-tmscore", str(args["min_tmscore"])]
        cmd += ["--output", str(out_dir)]

        result = _run(cmd, "Foldseek structural search")
        rj = _load_result_json(out_dir)

        if rj:
            top = rj.get("top_hit")
            summary = {
                "skill": "foldseek",
                "success": result["success"],
                "n_hits_raw": rj.get("n_hits_raw", 0),
                "n_hits_filtered": rj.get("n_hits_filtered", 0),
                "top_hit": top,
                "top_tmscore": top["tmscore"] if top else None,
                "novelty_verdict": (
                    "KNOWN (TM-score >= 0.9, prediction unnecessary)"
                    if top and top["tmscore"] >= 0.9
                    else "NOVEL (TM-score < 0.9, prediction recommended)"
                    if top
                    else "NO_HITS (no homologs found, structure is novel)"
                ),
                "top_hits": rj.get("hits", [])[:5],
                "output_dir": str(out_dir),
                "chat_summary": rj.get("chat_summary_lines", []),
            }
        else:
            summary = {
                "skill": "foldseek",
                "success": result["success"],
                "error": result.get("stderr", "")[:500] or result.get("error", ""),
                "stdout": result.get("stdout", "")[:500],
            }
        return json.dumps(summary, indent=2)

    # ── struct_predict ─────────────────────────────────────────────────────
    elif name == "struct_predict":
        out_dir = _skill_output_dir("struct_predictor")
        cmd = [sys.executable, str(SKILLS_DIR / "struct-predictor" / "struct_predictor.py")]
        input_yaml = args.get("input_yaml", "")
        if input_yaml:
            cmd += ["--input", input_yaml]
        else:
            cmd += ["--demo"]
        cmd += ["--output", str(out_dir)]

        result = _run(cmd, "Boltz-2 structure prediction")
        rj = _load_result_json(out_dir)
        report = _load_report_md(out_dir)

        summary = {
            "skill": "struct_predictor",
            "success": result["success"],
            "result_json": rj,
            "report_preview": report[:1000] if report else None,
            "output_dir": str(out_dir),
        }
        return json.dumps(summary, indent=2)

    # ── run_pharmgx ────────────────────────────────────────────────────────
    elif name == "run_pharmgx":
        out_dir = _skill_output_dir("pharmgx")
        input_file = args.get("input_file", "")
        demo = not bool(input_file)
        cmd = [sys.executable, str(CLAWBIO_PY), "run", "pharmgx"]
        if demo:
            cmd.append("--demo")
        else:
            cmd += ["--input", input_file]
        cmd += ["--output", str(out_dir)]

        result = _run(cmd, "PharmGx pharmacogenomics")
        report = _load_report_md(out_dir)
        return json.dumps({
            "skill": "pharmgx",
            "success": result["success"],
            "report_preview": report[:2000] if report else result.get("stdout", "")[:1000],
            "output_dir": str(out_dir),
        }, indent=2)

    # ── run_gwas_lookup ────────────────────────────────────────────────────
    elif name == "run_gwas_lookup":
        out_dir = _skill_output_dir("gwas")
        rsid = args.get("rsid", "")
        if not rsid:
            return json.dumps({"error": "rsid is required for run_gwas_lookup"})
        cmd = [sys.executable, str(CLAWBIO_PY), "run", "gwas",
               "--no-input-required", "--rsid", rsid, "--output", str(out_dir)]
        result = _run(cmd, f"GWAS lookup {rsid}")
        report = _load_report_md(out_dir)
        return json.dumps({
            "skill": "gwas_lookup",
            "rsid": rsid,
            "success": result["success"],
            "report_preview": report[:2000] if report else result.get("stdout", "")[:1000],
            "output_dir": str(out_dir),
        }, indent=2)

    # ── run_skill (generic) ────────────────────────────────────────────────
    elif name == "run_skill":
        skill = args.get("skill", "")
        if not skill:
            return json.dumps({"error": "skill name is required"})
        out_dir = _skill_output_dir(skill)
        input_file = args.get("input_file", "")
        extra_args = args.get("extra_args", [])
        cmd = [sys.executable, str(CLAWBIO_PY), "run", skill]
        if input_file:
            cmd += ["--input", input_file]
        else:
            cmd.append("--demo")
        cmd += ["--output", str(out_dir)]
        if extra_args:
            cmd += extra_args

        result = _run(cmd, f"ClawBio skill: {skill}")
        rj = _load_result_json(out_dir)
        report = _load_report_md(out_dir)
        return json.dumps({
            "skill": skill,
            "success": result["success"],
            "result_json": rj,
            "report_preview": report[:2000] if report else result.get("stdout", "")[:1000],
            "output_dir": str(out_dir),
        }, indent=2)

    else:
        return json.dumps({"error": f"Unknown tool: {name}"})


# ---------------------------------------------------------------------------
# Agent loop
# ---------------------------------------------------------------------------

MAX_TOOL_TURNS = 8  # prevent runaway loops


def run_agent(prompt: str, verbose: bool = True) -> str:
    """Run the agent loop for a single prompt. Returns the final text reply."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    messages: list[dict] = [
        {"role": "system", "content": SYSTEM_PROMPT},
        {"role": "user",   "content": prompt},
    ]

    if verbose:
        print(f"\n{'='*60}")
        print(f"Prompt: {prompt}")
        print(f"Model:  {MODEL}")
        print(f"Output: {OUTPUT_DIR}")
        print(f"{'='*60}")

    for turn in range(MAX_TOOL_TURNS):
        response = llm.chat.completions.create(
            model=MODEL,
            messages=messages,
            tools=TOOLS,
            tool_choice="auto",
        )

        msg = response.choices[0].message
        messages.append(msg)

        # No tool calls — LLM is done reasoning
        if not msg.tool_calls:
            final = msg.content or ""
            if verbose:
                print(f"\n{'─'*60}")
                print("Agent response:")
                print(final)
            return final

        # Execute each tool call
        for tc in msg.tool_calls:
            fn_name = tc.function.name
            try:
                fn_args = json.loads(tc.function.arguments)
            except json.JSONDecodeError:
                fn_args = {}

            if verbose:
                print(f"\n  Tool call: {fn_name}({json.dumps(fn_args, indent=2)})")

            tool_result = execute_tool(fn_name, fn_args)

            messages.append({
                "role": "tool",
                "tool_call_id": tc.id,
                "content": tool_result,
            })

    # Fallback if we hit MAX_TOOL_TURNS
    return "Agent reached maximum tool call limit. Check output directory for partial results."


# ---------------------------------------------------------------------------
# Interactive chat loop
# ---------------------------------------------------------------------------

def interactive_loop():
    """Simple REPL for multi-turn conversation."""
    global OUTPUT_DIR

    print("ClawBio Agent — interactive mode")
    print(f"Model: {MODEL}")
    print("Type your biological question. 'quit' to exit.\n")

    history: list[dict] = [{"role": "system", "content": SYSTEM_PROMPT}]

    while True:
        try:
            user_input = input("You: ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nBye.")
            break

        if user_input.lower() in ("quit", "exit", "q"):
            break
        if not user_input:
            continue

        # Fresh output dir per message
        OUTPUT_DIR = _ROOT / "output" / f"agent_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

        history.append({"role": "user", "content": user_input})
        messages = list(history)

        for turn in range(MAX_TOOL_TURNS):
            response = llm.chat.completions.create(
                model=MODEL,
                messages=messages,
                tools=TOOLS,
                tool_choice="auto",
            )
            msg = response.choices[0].message
            messages.append(msg)

            if not msg.tool_calls:
                reply = msg.content or ""
                print(f"\nAgent: {reply}\n")
                history.append({"role": "assistant", "content": reply})
                break

            for tc in msg.tool_calls:
                fn_name = tc.function.name
                try:
                    fn_args = json.loads(tc.function.arguments)
                except json.JSONDecodeError:
                    fn_args = {}
                print(f"  [calling {fn_name}...]")
                tool_result = execute_tool(fn_name, fn_args)
                messages.append({
                    "role": "tool",
                    "tool_call_id": tc.id,
                    "content": tool_result,
                })


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

_DEMO_PROMPT = (
    "You are my computational biologist. I want to understand if the Trp-cage "
    "miniprotein has any structurally similar proteins already known in the PDB. "
    "Run a structural homology search and tell me: is this fold novel, or does it "
    "match known structures? If it matches something, what is it and how similar is it?"
)


def main():
    global MODEL
    parser = argparse.ArgumentParser(
        description="ClawBio conversational agent — ask biological questions, get skill-backed answers"
    )
    parser.add_argument(
        "prompt",
        nargs="?",
        help="Biological question or analysis request",
    )
    parser.add_argument(
        "--interactive", "-i",
        action="store_true",
        help="Start interactive chat loop",
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Run with a canned demo prompt (structural homology question)",
    )
    parser.add_argument(
        "--model",
        default=MODEL,
        help=f"LLM model name (default: {MODEL})",
    )
    args = parser.parse_args()

    if args.model != MODEL:
        MODEL = args.model

    if args.interactive:
        interactive_loop()
    elif args.demo:
        run_agent(_DEMO_PROMPT)
    elif args.prompt:
        run_agent(args.prompt)
    else:
        parser.print_help()
        print("\nExample:")
        print('  python agent.py "Analyse my human genome and tell me if BRCA1 and TP53 share structural homology"')
        print("  python agent.py --demo")
        print("  python agent.py --interactive")


if __name__ == "__main__":
    main()
