"""
Launch a cavity MD run of the group code from presets in configs/.

The code is chosen with --code (see configs/codes.json): "cboamd" (the public
release, default) or "pdm" (polaritonic_deep_md, used for all runs before
2026-11-10). Each has its own clone next to this repo and its own conda env;
run.py uses that env's `cboamd`/`infrared`, whatever Python started it.

Builds in.json from  configs/defaults.json + configs/molecules/<mol>.json
+ configs/drivers/<driver>.json + command-line options, makes a fresh folder
under runs/<mol>/, records run_info.json (shared-repo commit, command, date),
then runs `cboamd` and `infrared` inside it. The code repos are never touched.

Examples:
    python scripts/run.py h2o pyscf                       # no cavity
    python scripts/run.py h2o pyscf --cavity --lam 0.1    # cavity on resonance
    python scripts/run.py co2 nep --cavity --lam 0.05 0.1 0.2   # coupling sweep
    python scripts/run.py co2 ff --steps 20 --dry-run     # just show in.json
    python scripts/run.py co2 pyscf --set basis=augccpvdz --name bigbasis
    python scripts/run.py co2 nep --cavity --lam 0.1 0.2 --campaign 2026-10-02_sweep
    python scripts/run.py h2o pyscf --cavity --lam 0.1 --code pdm   # old code
"""

import argparse
import datetime
import gzip
import json
import os
import platform
import shlex
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CONFIGS = REPO / "configs"
CM_PER_AU = 219474.63  # 1 hartree in cm^-1


def load(path):
    with open(path) as f:
        return {k: v for k, v in json.load(f).items() if not k.startswith("_")}


def choices(sub):
    return sorted(p.stem for p in (CONFIGS / sub).glob("*.json"))


def parse_value(text):
    """--set values: try JSON (numbers, lists, true/false), else keep as string."""
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return text


def code_info(name):
    """Repo folder, env and allowed drivers for a code in configs/codes.json."""
    codes = load(CONFIGS / "codes.json")
    info = {k: v for k, v in codes[name].items() if not k.startswith("_")}
    info["name"] = name
    info["repo"] = (REPO / info["repo"]).resolve()
    return info


def code_commit(code):
    def git(*args):
        return subprocess.run(["git", "-C", str(code["repo"]), *args],
                              capture_output=True, text=True).stdout.strip()
    return {"name": code["name"], "repo": str(code["repo"]), "env": code["env"],
            "commit": git("rev-parse", "--short", "HEAD"),
            "date": git("log", "-1", "--format=%ci"),
            "uncommitted_changes": bool(git("status", "--porcelain", "--untracked-files=no"))}


def tool(code, name):
    """Use the cboamd/infrared from the chosen code's env, not whatever is on PATH."""
    exe = Path(code["env"]) / "bin" / name
    if not exe.exists():
        sys.exit(f"'{name}' not found in {code['env']}/bin. Is {code['repo'].name} installed "
                 f"there? (pip install -e \"{code['repo']}\")")
    return str(exe)


def build_input(args, lam, code):
    mol = load(CONFIGS / "molecules" / f"{args.molecule}.json")
    jdata = load(CONFIGS / "defaults.json")
    jdata.update(load(CONFIGS / "drivers" / f"{args.driver}.json"))

    # Molecule-specific driver settings (e.g. NEP model files for CO2)
    jdata.update(mol.pop(args.driver, {}))
    for key in ("nep", "deepmd", "pyscf", "ff"):
        mol.pop(key, None)

    omega_cm = args.omega if args.omega is not None else mol.pop("omega_cm")
    mol.pop("omega_cm", None)
    jdata.update(mol)

    jdata["photons"] = args.cavity
    if args.cavity:
        jdata["omega_photon"] = [omega_cm / CM_PER_AU]
        jdata["lambda_photon"] = [lam]
    else:
        jdata.pop("lambda_vector", None)
    if args.steps is not None:
        jdata["steps"] = args.steps
    if args.dt is not None:
        jdata["timestep"] = args.dt
    for item in args.set:
        key, _, value = item.partition("=")
        jdata[key] = parse_value(value)

    # Paths in configs are relative to the repo; store them absolute so the run
    # folder works no matter where cboamd is started from.
    # "{code}" means the chosen code's repo (e.g. its bundled NEP models).
    for key in ("xyz_file", "nep_pot", "nep_dip", "nep_pol", "dp_pot", "dp_dip", "dp_pol"):
        if isinstance(jdata.get(key), str):
            jdata[key] = jdata[key].replace("{code}", str(code["repo"]))
        if isinstance(jdata.get(key), str) and not os.path.isabs(jdata[key]):
            jdata[key] = str(REPO / jdata[key])
            if not Path(jdata[key]).exists():
                print(f"warning: {key} not found: {jdata[key]}", file=sys.stderr)
    return jdata, omega_cm


def run_dir_name(args, lam, omega_cm):
    # Inside a campaign folder the date lives on the campaign, not each run
    parts = [args.driver] if args.campaign else [datetime.date.today().isoformat(), args.driver]
    if args.cavity:
        parts += ["cav", f"w{omega_cm:g}", f"lam{lam:g}"]
    if args.name:
        parts.append(args.name)
    return "_".join(parts)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("molecule", choices=choices("molecules"))
    ap.add_argument("driver", choices=choices("drivers"))
    ap.add_argument("--code", choices=[k for k in load(CONFIGS / "codes.json") if k != "default"],
                    default=load(CONFIGS / "codes.json")["default"],
                    help="which group code to run (configs/codes.json; default: %(default)s)")
    ap.add_argument("--cavity", action="store_true", help="couple to the cavity mode")
    ap.add_argument("--lam", type=float, nargs="+", default=[0.05],
                    help="coupling strength(s); several values = one run each")
    ap.add_argument("--omega", type=float, help="cavity frequency in cm^-1 (default: molecule preset)")
    ap.add_argument("--steps", type=int)
    ap.add_argument("--dt", type=float, help="timestep in fs")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                    help="override any in.json key, e.g. --set basis=augccpvdz")
    ap.add_argument("--name", help="extra tag appended to the run folder name")
    ap.add_argument("--campaign", help="group runs under runs/<molecule>/<campaign>/ (e.g. a sweep)")
    ap.add_argument("--no-ir", action="store_true", help="skip the infrared spectrum step")
    ap.add_argument("--dry-run", action="store_true", help="print in.json and folder, don't run")
    args = ap.parse_args()
    code = code_info(args.code)
    if args.driver not in code["drivers"]:
        sys.exit(f"driver '{args.driver}' is not available in --code {args.code} "
                 f"(has: {', '.join(code['drivers'])}).")

    lams = args.lam if args.cavity else [None]
    for lam in lams:
        jdata, omega_cm = build_input(args, lam, code)
        run_dir = REPO / "runs" / args.molecule / (args.campaign or "") / run_dir_name(args, lam, omega_cm)

        print(f"\n=== {run_dir.relative_to(REPO)}   [code: {args.code}]")
        if args.dry_run:
            print(json.dumps(jdata, indent=4))
            continue
        if run_dir.exists():
            sys.exit(f"{run_dir} already exists; use --name to make a new one.")
        run_dir.mkdir(parents=True)

        with open(run_dir / "in.json", "w") as f:
            json.dump(jdata, f, indent=4)
        with open(run_dir / "run_info.json", "w") as f:
            json.dump({"command": shlex.join([Path(sys.argv[0]).name, *sys.argv[1:]]),
                       "started": datetime.datetime.now().isoformat(timespec="seconds"),
                       "host": platform.node(),
                       "python": sys.executable,
                       "code": code_commit(code)}, f, indent=4)

        steps = [[tool(code, "cboamd"), "-i", "in.json"]]
        if not args.no_ir:
            steps.append([tool(code, "infrared"), "-m", "ase"])
        with open(run_dir / "run.log", "w") as log:
            for cmd in steps:
                print("$", " ".join(Path(c).name if i == 0 else c for i, c in enumerate(cmd)))
                proc = subprocess.Popen(cmd, cwd=run_dir, stdout=subprocess.PIPE,
                                        stderr=subprocess.STDOUT, text=True)
                for line in proc.stdout:
                    sys.stdout.write(line)
                    log.write(line)
                if proc.wait() != 0:
                    sys.exit(f"{cmd[0]} failed (exit {proc.returncode}); see {run_dir}/run.log")
        # Compressed copies are what git keeps (see .gitignore); the full outputs stay on disk
        for name in ("dipole.dat", "energy.dat"):
            with open(run_dir / name, "rb") as src, gzip.open(run_dir / f"{name}.gz", "wb") as dst:
                dst.write(src.read())
        print(f"done -> {run_dir.relative_to(REPO)}")
    if not args.dry_run:
        # keep the run -> code table (runs/CODES.md) current
        subprocess.run([sys.executable, str(REPO / "scripts" / "list_codes.py")])


if __name__ == "__main__":
    main()
