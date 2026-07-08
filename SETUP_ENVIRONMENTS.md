# SETUP_ENVIRONMENTS.md — Exact environment setup for all 3 computers

*Do these steps literally, in order. Commands are for **Windows PowerShell** (your main OS).
Where macOS/Linux differs, the line is marked `[mac/linux]`. Don't skip the verify steps —
they catch the failures that waste days later.*

Known facts about your repo (already true, 2026-07-07):
- GitHub remote already exists:
  `https://github.com/owen-hx-liu/Pan-Disease-Gene-Disease-Discovery-via-Evidence-Weighted-Knowledge-Graphs.git`
- `.gitignore` already excludes the big data (`data/`, `venv/`, `*.pkl`, embeddings) — good.
- The frozen graph `data/processed/edges_clean_integrated.csv` is **not** in git (correct);
  it must travel through a separate data channel (Step B below).
- Frozen graph hash: `8371ce0f...72cab5`, 360,549,016 bytes, 5,860,540 lines
  (see `ARTIFACT_HASHES.txt`).

---

## PART 0 — Who is who (assign this first)

Decide which physical computer is A, B, C and write it in `PROGRESS.md`.

| | Machine A | Machine B | Machine C |
|---|---|---|---|
| **Nickname** | GPU box | CPU box | Coordinator / laptop |
| **Must have** | NVIDIA RTX GPU (your 5070) | Lots of RAM (≥16 GB ideal) | Any laptop |
| **Installs** | Full stack **+ CUDA PyTorch** | Full stack (CPU) | Full stack (CPU) |
| **Runs** | KGE training, hybrid embeddings | baselines, leakage R0–R3, LOSO, null models | builds splits, figures, manuscript, Zenodo |
| **Git branch** | `machine-a` | `machine-b` | `machine-c` |

> If two machines have GPUs, both do Part 3. If only one does, that one is Machine A.

---

## PART A — One-time project prep (do ONCE, on the machine that has the repo now)

This is your current Windows machine. You're publishing the code + setting up the data channel
so B and C can join.

**A1. Make sure your latest code is on GitHub.**
```powershell
cd "C:\Users\owenh\Downloads\ScienceFairYear2"
git add -A
git commit -m "Add roadmap, scope, setup docs"
git push origin master
```
`[mac/linux]` same commands (no `.exe` differences here).

**A2. Create the shared DATA channel** (git can't hold the 360 MB graph). Pick ONE:
- **Google Drive (simplest):** make a folder named `SciencePaper-DATA`, right-click → Share →
  add the Google accounts of B and C as Editors.
- *(or)* a shared OneDrive / Dropbox / a USB SSD passed around / a home NAS.

**A3. Upload the frozen graph to that folder:**
upload `data\processed\edges_clean_integrated.csv` into `SciencePaper-DATA`.
(Later, Machine C also uploads the split files here.)

**A4. Confirm the hash is recorded** (already done for you in `ARTIFACT_HASHES.txt`). If you
ever regenerate the graph, recompute:
```powershell
Get-FileHash -Algorithm SHA256 data\processed\edges_clean_integrated.csv
```
`[mac/linux]` `sha256sum data/processed/edges_clean_integrated.csv`

✅ Part A done when: code is pushed to GitHub, and the graph file sits in a folder B and C can
download from.

---

## PART 1 — Base setup (do this on EVERY machine A, B, and C)

### Step 1 — Install Python 3.12 (exactly 3.12, not 3.13)
The pinned packages have 3.12 wheels. 3.13 may fail to install some of them.
```powershell
winget install -e --id Python.Python.3.12
```
- No winget? Download from https://www.python.org/downloads/release/python-3129/ and during
  install **tick "Add python.exe to PATH."**
- `[mac]` `brew install python@3.12`
- `[linux]` `sudo apt install python3.12 python3.12-venv python3.12-dev`

**Verify:**
```powershell
py -3.12 --version    # must print Python 3.12.x
```
`[mac/linux]` `python3.12 --version`

### Step 2 — Install Git
```powershell
winget install -e --id Git.Git
```
`[mac]` `brew install git`  ·  `[linux]` `sudo apt install git`

**Verify:** `git --version`

### Step 3 — Tell Git who you are (once per machine)
```powershell
git config --global user.name  "Your Name"
git config --global user.email "your-github-email@example.com"
```

### Step 4 — Clone the repo
```powershell
cd "$HOME\Documents"
git clone https://github.com/owen-hx-liu/Pan-Disease-Gene-Disease-Discovery-via-Evidence-Weighted-Knowledge-Graphs.git ScienceFairYear2
cd ScienceFairYear2
```
- First push/clone will ask you to log in to GitHub. If a browser popup doesn't appear, create
  a Personal Access Token (GitHub → Settings → Developer settings → Tokens → *classic*, scope
  `repo`) and paste it as the password.
- **Machine A only:** you already have the repo at `C:\Users\owenh\Downloads\ScienceFairYear2`
  — just `cd` there instead of cloning.

### Step 5 — Make and activate a virtual environment
```powershell
py -3.12 -m venv venv
.\venv\Scripts\Activate.ps1
```
- If PowerShell blocks activation with a script-policy error, run this once then re-activate:
  ```powershell
  Set-ExecutionPolicy -Scope Process -ExecutionPolicy RemoteSigned
  ```
- `[mac/linux]` `python3.12 -m venv venv && source venv/bin/activate`

**Verify:** your prompt now starts with `(venv)`.

### Step 6 — Upgrade pip
```powershell
python -m pip install --upgrade pip
```

### Step 7 — Get the frozen graph from the data channel
Download `edges_clean_integrated.csv` from the `SciencePaper-DATA` folder (Part A2) and place
it at exactly:
```
ScienceFairYear2\data\processed\edges_clean_integrated.csv
```
Create the folders if they don't exist:
```powershell
New-Item -ItemType Directory -Force data\processed | Out-Null
```
`[mac/linux]` `mkdir -p data/processed`  (then move the file in)

### Step 8 — VERIFY the graph hash (do NOT skip — this is the whole point of multi-machine work)
```powershell
Get-FileHash -Algorithm SHA256 data\processed\edges_clean_integrated.csv
```
`[mac]` `shasum -a 256 data/processed/edges_clean_integrated.csv`
`[linux]` `sha256sum data/processed/edges_clean_integrated.csv`

The output must equal the hash in `ARTIFACT_HASHES.txt`:
`8371ce0f2f837e68a7926ef86a3a3e2421b6ebadcc1d7f3b61dde17f0072cab5`
If it does **not** match → your download is corrupt/incomplete. Re-download. Do not proceed.

**Now go to your machine's specific part: A → Part 2, B → Part 3, C → Part 4.**

---

## PART 2 — Machine A only (the GPU box)

### A-1. Install / update the NVIDIA driver
- Install the **NVIDIA App** (or GeForce Experience) from https://www.nvidia.com/software/nvidia-app/
  and update to the latest Game-Ready/Studio driver (RTX 5070 = Blackwell, needs a recent driver).
- **Verify the GPU is visible:**
  ```powershell
  nvidia-smi
  ```
  You should see a table naming your RTX 5070 and a CUDA version ≥ 12.8. If `nvidia-smi` isn't
  found, the driver isn't installed correctly — fix before continuing.

### A-2. Install CUDA-matched PyTorch FIRST (before requirements.txt)
Your Blackwell GPU needs the CUDA 12.8 build:
```powershell
pip install torch==2.11.0 --index-url https://download.pytorch.org/whl/cu128
```

### A-3. Install the rest of the stack
```powershell
pip install -r requirements.txt
```
(torch is already satisfied from A-2, so pip won't replace it with the CPU build.)

### A-4. VERIFY the GPU is usable from Python
```powershell
python -c "import torch; print('CUDA OK:', torch.cuda.is_available(), '| GPU:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'NONE')"
```
Must print `CUDA OK: True | GPU: NVIDIA GeForce RTX 5070`. If it says `False`, the torch build
and driver don't match — reinstall torch with the cu128 index-url (A-2) and update the driver.

### A-5. VERIFY PyKEEN imports
```powershell
python -c "import pykeen; print('pykeen', pykeen.get_version())"
```
✅ Machine A done. Go to Part 5 (smoke test).

---

## PART 3 — Machine B only (the CPU box)

Machine B never trains KGE, so it does **not** need CUDA or a GPU. Two options:

**Option 1 (simplest — everything works):**
```powershell
pip install -r requirements.txt
```
This installs the CPU build of torch + pykeen too (a couple hundred MB you won't use, but no
import ever fails).

**Option 2 (leaner — skip the deep-learning packages B doesn't use):**
```powershell
pip install pandas==3.0.3 numpy==2.4.6 scipy==1.17.1 scikit-learn==1.9.0 networkx==3.6.1 matplotlib==3.10.9 tqdm==4.67.3
```
If any script later complains `ModuleNotFoundError: torch`, then run:
`pip install torch==2.11.0` (plain CPU build, no index-url).

**Verify the stack:**
```powershell
python -c "import pandas, numpy, scipy, sklearn, networkx, matplotlib, tqdm; print('CPU stack OK')"
```
✅ Machine B done. Go to Part 5 (smoke test).

---

## PART 4 — Machine C only (coordinator / laptop)

C builds the splits, makes figures, and writes — same CPU stack as B, plus release tooling.

### C-1. Install the CPU stack
```powershell
pip install -r requirements.txt
```
(or Option 2 from Part 3 to stay lean.)

### C-2. Make a free Zenodo account (for the final data/code DOI — needed in Phase 5)
- Go to https://zenodo.org → Sign up (you can log in with your GitHub account).
- Later you'll connect GitHub → Zenodo to snapshot a release. Nothing to install now; just
  have the account ready.

### C-3. Pick a manuscript tool (no install needed now, just decide)
- **Overleaf** (https://www.overleaf.com, free, LaTeX in the browser) — best for a real journal.
- or Google Docs / Word — fine for the JEI version.

**Verify the stack:**
```powershell
python -c "import pandas, numpy, scipy, sklearn, networkx, matplotlib, tqdm; print('coordinator stack OK')"
```
✅ Machine C done. Go to Part 5 (smoke test).

---

## PART 5 — Smoke test (run on ALL three machines before doing real work)

This proves the environment + the data are actually working end-to-end. Paste it as one block
in the activated `(venv)`:
```powershell
python -c "import pandas as pd, time; t=time.time(); df=pd.read_csv('data/processed/edges_clean_integrated.csv'); print('rows:', len(df)); print('cols:', list(df.columns)); print('relations:', df['relation'].nunique()); print('loaded in %.1fs' % (time.time()-t))"
```
Expected: `rows: 5860539`, columns `['source_id','relation','target_id','weight','dataset_sources']`,
`relations: 28`. If rows ≠ 5,860,539, your graph file is wrong (re-check Step 8 hash).

✅ If every machine prints the same `rows: 5860539`, all three environments are identical and
ready. You can now start the roadmap (Machine C builds `scripts/build_deleaked_splits.py` first).

---

## PART 6 — The daily git workflow (so 3 machines never overwrite each other)

**Each machine works on its own branch** (created once):
```powershell
git checkout -b machine-a      # machine A;  machine-b on B;  machine-c on C
git push -u origin machine-a
```

**Every time you start working:**
```powershell
git checkout master
git pull                       # get everyone's merged work
git checkout machine-a
git merge master               # bring your branch up to date
```

**Every time you finish a chunk (commit small, commit often):**
```powershell
git add -A
git commit -m "B: orthology blocker + R2 baseline eval"
git push
```

**To share your results with the others:** open a Pull Request on GitHub from your branch into
`master`, or just merge locally and push if you're comfortable. Small result files (JSON/CSV
tables, figures) go through git; **big files (embeddings, the graph) never do** — they go in the
`SciencePaper-DATA` channel with a hash added to `ARTIFACT_HASHES.txt`.

**The three rules that prevent disasters:**
1. **One owner per script.** A owns `kge_benchmark.py`; C owns `build_deleaked_splits.py`;
   B owns the leakage/LOSO scripts; C owns figures + manuscript. Don't edit someone else's file
   without telling them.
2. **Only Machine C regenerates the graph or the splits.** A and B consume them, never rebuild
   them — otherwise hashes drift and results stop being comparable.
3. **Verify the hash (Step 8) after every data-channel pull.** Always.

---

## Troubleshooting quick table

| Symptom | Fix |
|---|---|
| `Activate.ps1 cannot be loaded` | `Set-ExecutionPolicy -Scope Process RemoteSigned` then re-activate |
| `torch.cuda.is_available()` is False on A | reinstall torch with `--index-url .../cu128` (A-2); update NVIDIA driver |
| `nvidia-smi` not found | NVIDIA driver not installed / not on PATH |
| pip can't find a pinned version | you're on Python 3.13 — recreate the venv with `py -3.12` |
| hash mismatch at Step 8 | re-download the graph; the file is truncated/corrupt |
| `git push` asks for password endlessly | use a Personal Access Token (classic, `repo` scope) as the password |
| out-of-memory loading the graph on B/C | close other apps; ensure ≥8 GB free; use `dtype`/`usecols` in scripts |

---

*Once all three print `rows: 5860539` in Part 5, you're set up. Next: Machine C builds the
R0–R3 split generator (the backbone everything depends on) per `WORK_SPLIT_3_MACHINES.md`.*
