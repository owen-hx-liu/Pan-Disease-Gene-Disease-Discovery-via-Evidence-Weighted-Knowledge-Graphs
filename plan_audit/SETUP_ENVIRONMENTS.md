# SETUP_ENVIRONMENTS.md (audit plan) — environment setup for all 3 computers

*Windows PowerShell commands (your main OS). Lines marked `[mac/linux]` show the difference.
This is the old setup doc trimmed for the audit plan: same base, **+igraph**, **+Hetionet**,
and **KGE can run on free Colab/Kaggle** so the local GPU is no longer strictly required.
Don't skip verify steps.*

**What changed from the old setup doc**
- **Add `python-igraph`** to every machine (fast degree-preserving edge swaps for the R2 null —
  far faster than networkx).
- **Add a second graph to the data channel: Hetionet** (`hetionet-v1.0-edges.sif` + node file).
- **KGE (2 models only) may run on free Colab/Kaggle T4** instead of the local RTX — see Part 2.
- **Removed:** newer-release download tooling (that validation is cut in the audit plan).

Known facts about your repo (already true):
- GitHub remote:
  `https://github.com/owen-hx-liu/Pan-Disease-Gene-Disease-Discovery-via-Evidence-Weighted-Knowledge-Graphs.git`
- `.gitignore` excludes big data (`data/`, `venv/`, `*.pkl`, embeddings).
- Frozen Monarch graph `data/processed/edges_clean_integrated.csv`:
  hash `8371ce0f…72cab5`, 360,549,016 bytes, 5,860,540 lines (`ARTIFACT_HASHES.txt`).

---

## PART 0 — Who is who (assign first, write into PROGRESS.md)

| | Machine A | Machine B | Machine C |
|---|---|---|---|
| **Nickname** | GPU box | **CPU/RAM box (now the heavy one)** | Coordinator / laptop |
| **Must have** | RTX GPU (your 5070) *or* a Colab/Kaggle account | Lots of RAM (≥16 GB) + cores | Any laptop |
| **Installs** | Full stack + CUDA PyTorch *(or Colab)* | Full stack + **igraph** | Full stack + **igraph** |
| **Runs** | KGE (TransE, ComplEx) | **R2 degree-null, all baselines, Hetionet audit** | splits/regimes, figures, tables, manuscript, Zenodo |
| **Git branch** | `machine-a` | `machine-b` | `machine-c` |

> Role shift vs the old plan: KGE shrank to 2 small models on a filtered subgraph, so **A is no
> longer the bottleneck** — B is (the degree-permutation null over 5.86M edges). Keep B busy.
> If you have no usable local GPU at all, A becomes a second CPU box and KGE goes to Colab (Part 2B).

---

## PART A — One-time project prep (once, on the machine that has the repo)

**A1. Push code to GitHub.** (Same as before; first push may need `--force` once if the remote
still holds old raw-data uploads — see the original setup doc note. After that, normal `git push`.)
Never `git add` anything under `data/raw/`.

**A2. Create the shared DATA channel** (git can't hold the graphs). Pick ONE: Google Drive folder
`SciencePaper-DATA` shared to B and C as Editors, or OneDrive/Dropbox/NAS/USB-SSD.

**A3. Upload the two graphs to that folder:**
- `data\processed\edges_clean_integrated.csv` (Monarch, 360 MB) — primary.
- **Hetionet** (robustness graph) — get it once via **PART H below**, then drop both files into
  `SciencePaper-DATA`. ~2.25M edges / ~47k nodes — small.

**A4. Record hashes** in `ARTIFACT_HASHES.txt`:
```powershell
Get-FileHash -Algorithm SHA256 data\processed\edges_clean_integrated.csv
Get-FileHash -Algorithm SHA256 <path>\hetionet-v1.0-edges.sif
```
`[mac/linux]` `sha256sum <file>`

✅ Part A done when: code is on GitHub and **both graphs** sit in the data channel with hashes.

---

## PART H — Hetionet (the robustness graph) — get it ONCE

*Do this once on any machine with internet, then it travels through the data channel like the
Monarch graph. Hetionet is public, ~50 MB, human-only (no model organisms — that's why the
orthology regime R3 is Monarch-only and Hetionet runs R0–R2). Repo: https://github.com/hetio/hetionet*

### H1 — Download the two files
```powershell
New-Item -ItemType Directory -Force data\external\hetionet | Out-Null
$base = "https://github.com/hetio/hetionet/raw/main/hetnet/tsv"
Invoke-WebRequest "$base/hetionet-v1.0-nodes.tsv"    -OutFile data\external\hetionet\hetionet-v1.0-nodes.tsv
Invoke-WebRequest "$base/hetionet-v1.0-edges.sif.gz" -OutFile data\external\hetionet\hetionet-v1.0-edges.sif.gz
```
- If those URLs 404, swap `main` → `master` in `$base`.
- **Simplest alternative (branch-proof):** `git clone https://github.com/hetio/hetionet.git` — the
  files are then in `hetionet\hetnet\tsv\`; copy the two into `data\external\hetionet\`.
- `[mac/linux]` `curl -L -o hetionet-v1.0-edges.sif.gz "$base/hetionet-v1.0-edges.sif.gz"` (and the nodes file).

### H2 — Decompress the edges (Python is already installed; or skip — pandas reads `.gz` directly)
```powershell
python -c "import gzip,shutil; shutil.copyfileobj(gzip.open(r'data\external\hetionet\hetionet-v1.0-edges.sif.gz','rb'), open(r'data\external\hetionet\hetionet-v1.0-edges.sif','wb'))"
```
`[mac/linux]` `gunzip -k data/external/hetionet/hetionet-v1.0-edges.sif.gz`

### H3 — VERIFY it loaded and looks right
```powershell
python -c "import pandas as pd; df=pd.read_csv(r'data\external\hetionet\hetionet-v1.0-edges.sif', sep='\t'); print('edges:', len(df)); print('cols:', df.columns.tolist()); print('metaedges:', df['metaedge'].nunique()); print(df.head(3))"
```
Expected: **edges ≈ 2,250,197**, columns `['source','metaedge','target']`, ~24 distinct metaedges.
- The **prediction target for the audit = metaedge `DaG`** ("Disease–associates–Gene", ~12,623 edges).
- Genes are Entrez IDs (`Gene::1234`); diseases are DOID (`Disease::DOID:...`). Unit 9
  (`hetionet_audit.py`) maps `DaG` to the same `(gene, disease)` target format — it does **not**
  merge Hetionet with Monarch.

### H4 — Hash both files and append them to `ARTIFACT_HASHES.txt`
Paste this whole block (it computes the sha256 + byte count + line count for each file, prints the
lines, and appends them to the manifest in the exact `<sha256>  <bytes>  <lines>  <path>` format):
```powershell
$files = "data\external\hetionet\hetionet-v1.0-edges.sif", "data\external\hetionet\hetionet-v1.0-nodes.tsv"
Add-Content -Encoding utf8 ARTIFACT_HASHES.txt "`n# --- Hetionet v1.0 robustness graph (downloaded $(Get-Date -Format yyyy-MM-dd)) ---"
foreach ($f in $files) {
  $h = (Get-FileHash -Algorithm SHA256 $f).Hash.ToLower()
  $b = (Get-Item $f).Length
  $l = (Get-Content $f | Measure-Object -Line).Lines
  $line = "$h  $b  $l  $($f -replace '\\','/')"
  $line                                             # prints it so you can see it
  Add-Content -Encoding utf8 ARTIFACT_HASHES.txt $line
}
```
(The line count on the 2.25M-row edges file takes a few seconds — that's normal.) Then **upload the
two files to `SciencePaper-DATA`**, and commit the updated `ARTIFACT_HASHES.txt`. Every other machine
downloads them into `data\external\hetionet\` and re-verifies (Part 5). Only the machine that
downloaded Hetionet appends here — like the graph, hashes have one owner.

✅ Part H done when: both Hetionet files are hashed, in the data channel, and `H3` prints ~2.25M edges.

---

## PART 1 — Base setup (EVERY machine A, B, C)

### Step 1 — Python 3.12 (exactly 3.12)
```powershell
winget install -e --id Python.Python.3.12
py -3.12 --version        # must print 3.12.x
```
`[mac]` `brew install python@3.12` · `[linux]` `sudo apt install python3.12 python3.12-venv python3.12-dev`

### Step 2–3 — Git + identity
```powershell
winget install -e --id Git.Git
git config --global user.name  "Your Name"
git config --global user.email "your-github-email@example.com"
```

### Step 4 — Clone into a NON-OneDrive folder
> ⚠️ Not `Documents`/`Desktop` (OneDrive breaks venvs). Use `C:\dev` or `Downloads`.
```powershell
New-Item -ItemType Directory -Force C:\dev | Out-Null
cd C:\dev
git clone https://github.com/owen-hx-liu/Pan-Disease-Gene-Disease-Discovery-via-Evidence-Weighted-Knowledge-Graphs.git ScienceFairYear2
cd ScienceFairYear2
```
**Machine A** already has the repo at `C:\Users\owenh\Downloads\ScienceFairYear2` — just `cd` there.

### Step 5 — Virtual environment
```powershell
if (Test-Path venv) { Remove-Item -Recurse -Force venv }
py -3.12 -m venv venv
.\venv\Scripts\Activate.ps1     # if blocked: Set-ExecutionPolicy -Scope Process RemoteSigned, then re-activate
```
`[mac/linux]` `python3.12 -m venv venv && source venv/bin/activate`. Prompt must show `(venv)`.

### Step 6 — Upgrade pip
```powershell
python -m pip install --upgrade pip
```

### Step 7 — Get the graphs from the data channel
Place the Monarch graph at exactly `data\processed\edges_clean_integrated.csv`; place Hetionet at
`data\external\hetionet\` (both files).
```powershell
New-Item -ItemType Directory -Force data\processed | Out-Null
New-Item -ItemType Directory -Force data\external\hetionet | Out-Null
```

### Step 8 — VERIFY the Monarch graph hash (do NOT skip)
```powershell
Get-FileHash -Algorithm SHA256 data\processed\edges_clean_integrated.csv
```
Must equal `8371ce0f2f837e68a7926ef86a3a3e2421b6ebadcc1d7f3b61dde17f0072cab5`
(from `ARTIFACT_HASHES.txt`). Mismatch → re-download; do not proceed. Verify Hetionet's hash too.

**Then go to your machine's part: A → Part 2, B → Part 3, C → Part 4.**

---

## PART 2 — Machine A (KGE)

You have two options. **2A (local RTX)** is faster if the driver cooperates; **2B (Colab/Kaggle)**
is the safety net and is genuinely fine here because the models are tiny (dim 64, filtered subgraph).

### 2A — Local RTX 5070
```powershell
nvidia-smi                                   # must list the 5070, CUDA ≥ 12.8
pip install torch==2.11.0 --index-url https://download.pytorch.org/whl/cu128
pip install -r requirements.txt
pip install python-igraph                     # for any local graph ops
python -c "import torch; print('CUDA OK:', torch.cuda.is_available())"   # must be True
python -c "import pykeen; print('pykeen', pykeen.get_version())"
```
If `CUDA OK` is False, reinstall torch with the cu128 index-url and update the NVIDIA driver.

### 2B — Free cloud GPU (Colab or Kaggle) — use if no local GPU or the driver fights you
Nothing to install locally. In a Colab/Kaggle notebook with a **T4** runtime:
```python
!pip -q install pykeen==1.11.1 torch --upgrade
# upload data/processed/splits/*.csv (or mount Drive), then run scripts/run_kge.py
```
Train on the **filtered subgraph** (see BUILD_GUIDE Unit K) so each run is tens of minutes.
Download the resulting `data/processed/results/kge_*.json` back into the repo and commit.

✅ Machine A done → Part 5.

---

## PART 3 — Machine B (CPU/RAM — the heavy machine now)

B runs all baselines, the R2 degree-permutation null, and the Hetionet audit. No GPU needed.
```powershell
pip install pandas numpy scipy scikit-learn networkx matplotlib tqdm python-igraph
python -c "import pandas,numpy,scipy,sklearn,networkx,igraph,matplotlib,tqdm; print('CPU+igraph OK')"
```
(or `pip install -r requirements.txt` then `pip install python-igraph`.) **igraph is mandatory on B** —
the degree-preserving double-edge swaps over 5.86M edges are impractically slow in networkx.

✅ Machine B done → Part 5.

---

## PART 4 — Machine C (coordinator / laptop)

C builds the splits/regimes, makes figures/tables, writes, and releases.
```powershell
pip install pandas numpy scipy scikit-learn networkx matplotlib tqdm python-igraph
python -c "import pandas,numpy,scipy,sklearn,networkx,igraph,matplotlib,tqdm; print('coordinator OK')"
```
- **Zenodo:** make a free account (log in with GitHub) — needed for the final data/code DOI.
- **Manuscript tool:** Overleaf (LaTeX, best for a real journal) or Google Docs.

✅ Machine C done → Part 5.

---

## PART 5 — Smoke test (ALL three machines, in the activated venv)

```powershell
python -c "import pandas as pd, time; t=time.time(); df=pd.read_csv('data/processed/edges_clean_integrated.csv'); print('rows:', len(df)); print('cols:', list(df.columns)); print('relations:', df['relation'].nunique()); print('loaded in %.1fs' % (time.time()-t))"
```
Expected `rows: 5860539`, cols `['source_id','relation','target_id','weight','dataset_sources']`,
`relations: 28`. If rows differ, the graph file is wrong (re-check Step 8).

✅ When all three print `rows: 5860539` and import igraph, environments are identical. Start the
build plan (Machine C builds the R0–R3 regime generator first — BUILD_GUIDE Unit 2).

---

## PART 6 — Daily git workflow (unchanged)

Each machine on its own branch (`machine-a/-b/-c`); merge to `main` via PRs. Start each session
`git checkout main; git pull; git checkout <mine>; git merge main`. Commit small/often; push.
Small result files (JSON/CSV/figures) go in git; **big files (graphs, embeddings) never do** — data
channel + a hash in `ARTIFACT_HASHES.txt`.

**Three rules that prevent disasters:**
1. **One owner per script.**
2. **Only Machine C regenerates the graph or the splits/regimes.** A and B consume, never rebuild.
3. **Verify the hash after every data-channel pull.** Always.

---

## Troubleshooting (delta from the old table)

| Symptom | Fix |
|---|---|
| `import igraph` fails | `pip install python-igraph` (note: package is `python-igraph`, module is `igraph`). |
| igraph edge-swap is slow / OOM on B | do the swap on the *undirected simple* graph, chunk swaps, keep int32 ids; see BUILD_GUIDE Unit R2. |
| No local CUDA on A | use Part 2B (Colab/Kaggle) — the audit plan is designed to not need the local GPU. |
| Hetionet columns don't match | Hetionet `.sif` is `source⇥metaedge⇥target`; map its Gene/Disease metaedges in Unit H, don't merge it with Monarch. |
