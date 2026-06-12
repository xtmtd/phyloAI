# Checkpoint and Resume Design Specification

**Date:** 2026-06-12  
**Status:** Draft for user review  
**Parent spec:** `2026-06-07-phyloai-design.md`

---

## 1. Purpose

Long-running PhyloAI commands can run for hours or days. A run may stop because of power loss, scheduler timeout, external-tool failure, machine reboot, or user interruption. The original design deferred resume support, which leaves expensive partial analyses hard to recover.

This specification adds a shared checkpoint and explicit `--resume` model for long-running commands. The first implementation target is `phyloai pretree align`; later commands should reuse the same core mechanism.

The design goals are:

- preserve the existing safe default: non-empty output directories still fail unless the user explicitly chooses `--overwrite` or `--resume`
- make resume reproducible: resumed runs must use exactly the same parameters as the original run
- avoid duplicating `result.json`: checkpoints store only the minimum state needed to continue work
- support both batch-per-file commands and tool-native resume commands such as IQ-TREE and PhyloBayes

---

## 2. Scope

### 2.1 In scope

- Shared checkpoint file format: `checkpoint.json`
- Shared checkpoint helpers in `phyloai/core/checkpoint.py`
- New `--resume` flag for long-running pipeline commands
- First command adoption: `phyloai pretree align`
- Main design update to replace the old "resume later" note in Section 9.5

### 2.2 Out of scope

- Background job daemon or persistent task queue
- Cross-command DAG orchestration for `phyloai run`
- Automatic resume without explicit `--resume`
- Resume support for short utility commands (`pretree convert`, `pretree stats`) in the first implementation
- Partial parameter compatibility rules; all parameters must match exactly for now

---

## 3. User-Facing CLI Semantics

Long-running commands support three mutually exclusive output-directory modes:

| Mode | Behavior |
|------|----------|
| default | If `--output-dir` exists and is non-empty, exit 1 with a clear message |
| `--overwrite` | Delete and recreate `--output-dir`, then start from scratch |
| `--resume` | Load `checkpoint.json`, validate parameters, skip verified successes, rerun incomplete or failed tasks |

`--overwrite` and `--resume` are mutually exclusive. Supplying both exits with code 1.

`--resume` requires `checkpoint.json` inside the output directory. If the file is missing, malformed, or incompatible with the command, the command exits with code 1 and tells the user to use `--overwrite` if they want a clean restart.

`--resume` requires the current resolved parameters to match the checkpoint parameters exactly. This includes analysis parameters and run-control parameters such as `--threads` and `--quiet`. Exact matching is intentionally strict in the first version to avoid silently mixing outputs from different analytical settings.

If the checkpoint indicates that the command already completed successfully, `--resume` does not rerun work. The command reports that the run is already complete and points to the existing `result.json`.

With `--dry-run --resume`, the command does not execute tools. It reports how many tasks would be skipped, rerun, or considered invalid after output verification.

---

## 4. Checkpoint File Model

Each resumable command writes `checkpoint.json` in its output directory, alongside `result.json` and the command log.

Example layout for `pretree align`:

```text
runs/pretree/align/
├── seqs/
├── align.log
├── checkpoint.json
└── result.json
```

`result.json` remains the final structured command result. `checkpoint.json` is a recovery manifest, not a reporting artifact.

### 4.1 Top-level schema

```json
{
  "schema_version": 1,
  "step": "pretree.align",
  "command": "phyloai pretree align --seq-dir raw --method linsi ...",
  "status": "running",
  "params_hash": "sha256:...",
  "params": {},
  "started_at": "2026-06-12T10:00:00",
  "updated_at": "2026-06-12T10:30:00",
  "completed_at": null,
  "tasks": []
}
```

Top-level `status` values:

| Status | Meaning |
|--------|---------|
| `running` | Command started and may have incomplete tasks |
| `success` | Command completed and wrote final `result.json` |
| `error` | Command stopped after one or more unrecoverable command-level errors |
| `interrupted` | Command stopped before normal completion, for example by signal or scheduler timeout |

### 4.2 Task schema

Tasks store only the state needed to decide whether to skip or rerun that unit. They do not duplicate per-task metrics that belong in `result.json`.

```json
{
  "task_id": "gene1",
  "status": "success",
  "input": "raw/gene1.fa",
  "outputs": {
    "aa": "runs/pretree/align/seqs/gene1.fa",
    "nt": null
  },
  "attempts": 1,
  "reason": null,
  "updated_at": "2026-06-12T10:30:00"
}
```

Task `status` values:

| Status | Resume behavior |
|--------|-----------------|
| `pending` | Run task |
| `running` | Treat as interrupted and rerun task |
| `success` | Verify outputs; skip if valid, rerun if invalid |
| `failed` | Rerun task by default |
| `skipped` | Keep skipped unless the command's normal scanning or validation now classifies it differently |

The `outputs` object is command-specific but should contain only paths needed for resume verification. Detailed metrics such as `n_taxa`, `alignment_length`, `wall_time`, warnings, likelihoods, and scores belong in final outputs or command-specific partial result files, not in the shared checkpoint.

---

## 5. Atomic Writes and Crash Safety

Checkpoint writes must be atomic:

1. Serialize the full checkpoint to `checkpoint.json.tmp`
2. Flush the file
3. Replace `checkpoint.json` with the temp file using an atomic filesystem replace (`os.replace`)

`os.replace` is atomic on POSIX and Windows, so a reader never sees a half-written
`checkpoint.json` even if the process is killed mid-write. This protects against
`ctrl+C`, `SIGTERM`, and process crashes **without** an `fsync`.

### 5.1 When to call `fsync`

`fsync` only adds durability against **hard power loss / kernel panic** (it forces
the bytes from the OS page cache onto the physical disk). It does **not** improve
safety against process termination, which `os.replace` already handles.

`fsync` on a large file is expensive. Measured cost on the reference machine for a
~450 KB checkpoint (1066 tasks) was **~258 ms per write**. Calling it after every
task is therefore **O(N) writes × O(N) file size = O(N²)** of pure I/O — on the
order of **~275 s of overhead** for ~1000 genes, which made resumed and even fresh
runs dramatically slower than expected. See Section 14.2.

Therefore:

- During the task loop, write **without** `fsync` (atomicity comes from `os.replace`).
- Call `fsync` **only** on the final completion write and on the interrupt-handler write.

`save_checkpoint_atomic(checkpoint, path, *, fsync: bool = False)` exposes this as a
parameter; default is `False`.

### 5.2 Write frequency: throttle, do not write per task

Writing the full checkpoint after every single task is both the O(N²) problem above
and unnecessary for correctness. Instead:

- Flush at most once every `CHECKPOINT_FLUSH_INTERVAL` seconds (currently `2.0`)
  during the loop, without `fsync`.
- On `KeyboardInterrupt`, force one final flush, set top-level `status = "interrupted"`,
  then re-raise.
- On normal completion, set `status = "success"` and do one `fsync` flush.

The worst case lost progress on interrupt is bounded by the flush interval (a few
seconds of finished tasks), which resume will simply rerun. This is an acceptable
trade for removing per-task `fsync` cost.

Commands logically update the in-memory checkpoint at these points; the *disk flush*
is throttled as above:

- after initial task discovery (flush once, this is cheap and rare)
- when a task moves to `running`
- when a task reaches `success`, `failed`, or `skipped`
- when command-level status changes to `success`, `error`, or `interrupted`
  (always flush immediately for these terminal transitions)

For highly parallel commands, per-task updates are applied in the parent process as
workers complete. Worker processes must never write the checkpoint directly.

---

## 6. Parameter Matching

Each resumable command builds a resolved parameter dictionary after applying defaults and config-file values. The checkpoint stores this dictionary and a stable SHA-256 hash of its canonical JSON representation.

On resume:

1. Rebuild the resolved parameter dictionary from the current invocation
2. Compute the hash using the same canonicalization rules
3. Compare it with `params_hash`
4. If different, exit 1 and show a concise mismatch message

The first version requires exact matching. If users want to change any parameter, they must start over with `--overwrite`.

---

## 7. `pretree align` Resume Flow

### 7.1 Initial run

On a non-resume run, `pretree align`:

1. Applies the normal output conflict policy
2. Scans `--seq-dir`
3. Resolves `--seq-type auto`, tool paths, and other parameters
4. Creates `checkpoint.json` with one task per valid input gene
5. Runs alignment tasks in parallel
6. Updates each task status after completion
7. Builds final `result.json` from current output files and task outcomes
8. Marks the checkpoint top-level status as `success` after `result.json` is written

### 7.2 Resume run

On `--resume`, `pretree align`:

1. Loads `checkpoint.json`
2. Verifies step name and schema version
3. Verifies exact parameter hash match
4. Validates outputs for tasks marked `success`
5. Skips successful tasks whose outputs are still valid
6. Reruns tasks with status `pending`, `running`, or `failed`
7. Reruns tasks whose recorded success outputs are missing or invalid
8. Rebuilds final `result.json` from all valid outputs and current task statuses

AA output validation uses the shared sequence-output validation helper with `require_aligned=True`. In backtranslation mode, recorded NT outputs are also checked for existence, parseability, and aligned sequence lengths.

An AA success with missing or invalid NT output is treated as incomplete for that task. The first implementation reruns the whole gene task for simplicity and correctness. A later optimization may rerun only the backtranslation substep if whole-task reruns become too expensive.

### 7.4 A task must mean "all of its declared outputs exist" — atomic task units

**Critical invariant:** a task may only be marked `success` once **every output it
declares** in the checkpoint (`outputs.aa`, `outputs.nt`, …) has been produced and
will pass the resume verifier. If a task declares an `nt` output but is marked
`success` before the NT file exists, resume verification will fail and rerun it — so
the checkpoint is actively lying about what is done.

This drove the original `--backtrans` resume bug (Section 14.1): AA alignment and NT
backtranslation were two **separate phases** (align all genes in the pool, *then*
backtranslate all genes), but each gene's checkpoint task declared both `aa` and `nt`
outputs. Genes got marked `success` after the AA phase, before any NT existed.
Interrupting during the (long) AA phase left a checkpoint full of `success` tasks
whose `nt` files did not exist, so resume reran everything.

**Design rule:** the unit of checkpointing must equal the unit of completed work.
If a "gene" task owns both an AA and an NT output, then all stages that produce those
outputs must complete **for that gene** before it is marked `success`. Concretely,
run the gene's backtranslation **inline, immediately after its AA alignment
completes**, within the same worker-completion handler — not in a separate later
phase. Each gene is then a single atomic, resumable unit.

If you genuinely need multi-phase processing where phases cannot be fused, then either:

1. split the work into separate checkpoint tasks per phase (e.g. `gene1.aa`,
   `gene1.nt`) so each phase is independently tracked, or
2. add an intermediate task status (e.g. `aa_done`) and make the resume verifier and
   planner phase-aware.

Do **not** declare an output in a task and then mark that task `success` before the
output exists.

### 7.3 Final result reconstruction

The final `result.json` is reconstructed from:

- checkpoint task states
- validated output files
- normal command-level scan results
- command-specific output readers that compute fields such as `n_taxa` and `alignment_length`

This keeps checkpoint files small and prevents stale metrics from surviving after output files change.

---

## 8. Future Command Mapping

The checkpoint framework should be command-agnostic. Each command defines its task granularity and output-verification function.

| Command | Task unit | Resume strategy |
|---------|-----------|-----------------|
| `pretree trim` | one MSA file | skip verified trimmed alignments; rerun failed/incomplete files |
| `pretree metrics` | one MSA or MSA/tree pair | skip metric records that exist in command-specific partial outputs; rerun missing records |
| `tree genetree` | one gene alignment | skip verified tree files; rerun failed/incomplete gene trees |
| `tree iqtree` | one matrix analysis | prefer IQ-TREE native checkpoint/resume; PhyloAI records tool state and verifies final outputs |
| `tree phylobayes` | one chain | use PhyloBayes chain resume; PhyloAI tracks chain status and convergence-check outputs |
| `posttree concordance` | tree/branch or gene-tree batch unit | skip verified output tables; rerun incomplete units |
| `posttree topology` | one hypothesis/test unit | skip verified test outputs; rerun incomplete tests |
| `posttree dating` | one MCMCTree run or chain | use tool-native restart where available; verify final posterior outputs |
| `posttree signal` | one site/gene/hypothesis unit | skip verified support tables; rerun incomplete units |
| `posttree syserror` | one atomic diagnostic operation | checkpoint each atomic operation result |
| `posttree simulate` | one replicate or replicate batch | skip verified simulation outputs; rerun missing/failed replicates |

For commands backed by tools with native resume, PhyloAI should not replace the tool's checkpointing. It should wrap it: record the command state, call the tool with the correct resume flags or existing working directory, and verify final outputs.

---

## 9. Core API Design

Create `phyloai/core/checkpoint.py` with small, reusable helpers.

Proposed data model:

```python
@dataclass
class CheckpointTask:
    task_id: str
    status: str
    input: str
    outputs: dict[str, str | None]
    attempts: int = 0
    reason: str | None = None
    updated_at: str | None = None


@dataclass
class Checkpoint:
    schema_version: int
    step: str
    command: str
    status: str
    params_hash: str
    params: dict[str, Any]
    started_at: str
    updated_at: str
    completed_at: str | None
    tasks: list[CheckpointTask]
```

Proposed helper functions:

```python
def canonical_params_hash(params: dict[str, Any]) -> str: ...
def load_checkpoint(path: Path) -> Checkpoint: ...
def save_checkpoint_atomic(checkpoint: Checkpoint, path: Path, *, fsync: bool = False) -> None: ...
def validate_resume_params(checkpoint: Checkpoint, params: dict[str, Any]) -> None: ...
def summarize_resume_tasks(checkpoint: Checkpoint, verifier: Callable[[CheckpointTask], bool]) -> dict[str, int]: ...
```

`save_checkpoint_atomic` takes `fsync` (default `False`) per Section 5.1: pass
`fsync=True` only on terminal/interrupt writes.

The core module should not know command-specific output formats. Commands provide verifier functions and result reconstruction logic.

---

## 10. Error Handling

Resume-specific failures use exit code 1 because they are user/input state errors:

- `--resume` with no checkpoint
- malformed checkpoint JSON
- checkpoint step does not match the current command
- unsupported checkpoint schema version
- parameter hash mismatch
- `--resume` combined with `--overwrite`

External tool failures keep the existing command semantics. For batch commands, at least one successful task can still produce exit 0 with failed tasks recorded, while all tasks failing exits non-zero according to the command's normal failure policy.

If a process receives an interrupt that the CLI can catch, it should mark the checkpoint top-level status as `interrupted` before exiting. Hard power loss will leave the top-level status as `running`; resume treats old `running` tasks as interrupted and reruns them.

---

## 11. Documentation Requirements

Implementation must update:

- `docs/commands/pretree-align.md` to document `--resume`, checkpoint behavior, and parameter mismatch errors
- `README.md` command index or common CLI behavior section if one exists
- main design Section 9.5 to replace "`--resume` is not in scope" with the explicit checkpoint policy

Future command docs should include a short "Resume behavior" section whenever they support `--resume`.

---

## 12. Testing Requirements

Shared checkpoint tests:

- atomic save writes valid JSON
- load rejects malformed JSON with a clear error
- parameter hash is stable for equivalent dictionaries
- parameter mismatch is detected
- unsupported schema version is rejected

`pretree align` tests:

- default non-empty output directory still errors
- `--overwrite` and `--resume` are mutually exclusive
- `--resume` without checkpoint errors
- parameter mismatch errors
- successful tasks with valid outputs are skipped
- failed tasks are retried by default
- `running` tasks from an interrupted run are retried
- success tasks with missing or invalid outputs are retried
- final `result.json` includes both previously completed and newly completed tasks
- dry-run resume reports skip/rerun counts without executing tools
- **multi-output task interrupt regression (Section 14.1):** for a command with
  multiple declared outputs per task (e.g. `--backtrans` producing AA + NT), an
  interrupted run must NOT re-process tasks whose outputs already exist. The test
  must run real tools, complete a run, simulate an interrupt by reverting some tasks
  to `running` and deleting their outputs, resume, and assert that completed tasks are
  **not** re-executed (e.g. by checking output file mtimes are unchanged). A unit test
  that only marks tasks `success` without exercising the real two-stage execution path
  will not catch this class of bug.

---

## 13. Main Design Updates

The parent spec should be updated as follows:

- Section 9.2 shared parameter registry: add `--resume` for long-running pipeline commands
- Section 9.5 output directory conflict policy: replace the old statement that resume is out of scope with the explicit default/overwrite/resume behavior from this spec
- Section 9.6 logging: note that resumed runs append to the existing command log and should include resume context in new log entries

This spec does not change the development phase order. It adds shared infrastructure that should be implemented before future long-running subcommands are added, with `pretree align` serving as the first adoption point.

---

## 14. Lessons Learned (post-implementation)

This section records two real bugs found while testing the first `pretree align`
implementation, so future adopters avoid repeating them. Both produced the same
user-visible symptom — **`--resume` re-ran everything and was no faster (or slower)
than a fresh run** — but had different root causes.

### 14.1 "Success" was recorded before all declared outputs existed

**Symptom:** Interrupt a `--backtrans` run partway through, then `--resume`. All
already-aligned genes were re-aligned from scratch; wall time matched a full run.

**Root cause:** Execution was two-phase — the parallel pool aligned *all* genes (the
expensive part), and a separate later loop ran backtranslation for *all* genes. But
each gene's checkpoint task declared both an `aa` and an `nt` output, and genes were
marked `success` at the end of the AA phase. Interrupting during the AA phase left a
checkpoint where many tasks were `success` but their `nt` files did not exist. On
resume, the output verifier (which checks NT too) returned `False` for every such
task, so the planner requeued them all.

Evidence that pinpointed it: in the post-resume checkpoint, completed-then-rerun
tasks had `attempts == 4` (run1 running+success, run2 running+success) while
interrupted-mid-AA tasks had `attempts == 3`. The `attempts == 4` group proved tasks
recorded `success` in run1 were still rerun in run2 — i.e. verification was failing on
existing-but-incomplete output sets.

**Fix:** Fuse the stages. Run each gene's backtranslation inline, immediately after
its AA alignment completes, in the same worker-completion handler, and only then mark
the gene `success`. Each gene is now one atomic resumable unit. See Section 7.4.

**Generalizable rule:** *The checkpoint task is a promise that all its declared
outputs exist and are valid. Never mark `success` until that promise is true.* If
work is inherently multi-phase, either split it into per-phase tasks or make the
verifier/planner phase-aware — do not record premature success.

### 14.2 Per-task `fsync` made checkpointing O(N²)

**Symptom:** After fixing 14.1, resume worked correctly but the whole run (fresh and
resumed) became far slower than the original non-checkpointed implementation.

**Root cause:** The checkpoint was saved after every task with an `fsync`. The file
grows with the number of tasks, so writing it N times, each O(file size), is O(N²).
For ~1066 genes the checkpoint was ~450 KB and each `fsync` write measured ~258 ms —
roughly **275 s** of pure checkpoint I/O added to every run.

**Fix:** Two changes (Section 5):

1. Drop `fsync` from the hot path. `os.replace` already guarantees atomicity against
   process termination; `fsync` only protects against power loss and is the expensive
   part. Keep `fsync` only on the final completion write and the interrupt write.
2. Throttle disk flushes to at most once every `CHECKPOINT_FLUSH_INTERVAL` (2 s)
   during the loop, forcing a flush on interrupt and at completion.

After this, per-task overhead dropped from ~258 ms to effectively zero on the hot
path; worst-case lost progress on interrupt is a couple of seconds of finished tasks,
which resume cheaply reruns.

**Generalizable rules:**

- Distinguish *atomicity* (use `os.replace`, cheap) from *durability* (`fsync`,
  expensive). For interrupt/`ctrl+C` safety you only need atomicity.
- Beware O(N²) when you serialize a growing whole-file state once per unit of work.
  Throttle by time (or size delta), and only force a durable write at terminal points.

### 14.3 Debugging method that worked

The bugs were found by inspecting the **actual on-disk `checkpoint.json` after a real
interrupt** (task status counts, `attempts` distribution, and whether declared output
files existed), and by **measuring** the atomic-write cost directly rather than
guessing. Unit tests that only set task status to `success` in memory passed the whole
time because they never exercised the real two-phase execution or the real file I/O.
A faithful regression test must drive the real execution path and simulate a true
interrupt (see Section 12).
