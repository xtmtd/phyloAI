# PhyloAI `tree bi pb` Stop-Signal Repair Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: use `superpowers:subagent-driven-development` or `superpowers:executing-plans` to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

> **Commits are user-driven only.** No task runs `git add` or `git commit`. Each task ends with “Report done — no commit”. The user reviews the diff and commits when they ask. A complete patch is exported after each test-passing stage.

**Goal:** Replace PR #1’s provisional `result.json` shutdown path with a CLI-owned `StopController`, a library-owned `stop_state.json` sidecar, atomic result publication, deterministic soft/forced stop semantics, and matching tests/docs.

**Baseline:** PR #1 head `36d4f96`, base `ddbfb98`. Existing BI + MCP + CLI tests (five files): **85 collected**. The six-file baseline command adds `tests/core/test_schema.py` (11 tests); expected **96 passed** before this change. `tests/tree/test_bi_readpb.py` has **6 pre-existing failures** on both PR and base; do not fix them in this plan.

**Spec references:**

- `docs/superpowers/specs/2026-06-07-phyloai-design.md` §9.3
- `docs/superpowers/specs/2026-06-23-phyloai-tree-bi-design.md` §3.1, §3.4, §5.2, §6, §7, §9, §10
- `docs/superpowers/specs/2026-06-21-phyloai-json-output-standard.md` §4, §7
- `docs/superpowers/specs/2026-06-27-phyloai-ai-integration-design.md` §2.3
- `docs/superpowers/specs/2026-07-26-phyloai-tree-bi-subcommands-design.md` and its plan’s “CLI sole writer” convention
- `docs/superpowers/specs/2026-06-12-checkpoint-resume-design.md` §5

**Architecture decisions already fixed:**

- CLI owns the single `StopController`; library never installs/restores handlers on the CLI path.
- Library writes only `stop_state.json`; CLI remains sole writer of `result.json`.
- First stop: soft stop → final success/0 with `shutdown_state: stopped`.
- Second stop before commit: forced stop → error + `128+signal`, `shutdown_state: forced`. Cleanup targets the direct chain process and is best-effort; real MPI descendants depend on `mpirun`/scheduler process-group handling.
- `begin_commit()` sets `phase = committing` first, then checks `pending_request = stop_level > handled_level`; historical requests do not block commit.
- Only the CLI final publisher calls `begin_commit()`. A late pending request after library return is handled by the library finalization API using the existing run context; never rerun analysis.
- Commit failure (backup/serialize/flush/fsync/replace) changes an original success 0 to 4; an already-selected forced code 130/143/129 is preserved; existing 1/2/3 handling is unchanged.
- Direct library calls use `StopController(install_handlers=False)` plus a unified-loop `KeyboardInterrupt` boundary; first Ctrl+C soft, second forced.
- Failed-chain cleanup uses the same bounded wait/controller checks; the original chain-failure message is preserved.
- Stale `stop_state.json` is cleared only after input preparation (including matrix conversion) and chain-command validation succeed.
- `check_status` gives an existing sidecar priority; a corrupt sidecar returns `unknown` + `stop_state_error`, never a stale success.
- Version bump `0.5.1` is a separate change.

---

## File Map

| File | Action | Responsibility |
|---|---|---|
| `phyloai/tree/bi.py` | Modify | `StopController`, `stop_state.json`, soft/forced shutdown, failed-chain bounded cleanup, late finalization API, direct-call fallback |
| `phyloai/cli/commands/tree.py` | Modify | Controller lifecycle, backup timing, atomic final/error publish, business-code preservation, exit 4, sidecar cleanup |
| `phyloai/core/schema.py` | Modify | Shared atomic JSON helper, permission preservation |
| `phyloai/mcp/tools/utils.py` | Modify | Sidecar priority, corrupt-sidecar `unknown`, tolerant reads |
| `tests/tree/test_bi.py` | Modify | Controller/SIG_IGN, direct double Ctrl+C, failed-chain force, sidecar timing |
| `tests/cli/test_tree_bi.py` | Modify | Backup timing, TypeError/ValueError publication, business-code preservation, error-write masking |
| `tests/cli/test_tree_bi_signals.py` | Modify | SIGTERM/SIGHUP/nohup, sidecar hard kill, forced exit codes, normal completion |
| `tests/core/test_schema.py` | Modify | Atomic write, old-target preservation, permissions |
| `tests/mcp/test_utils.py` | Modify | Sidecar priority, corrupt sidecar, running sidecar |
| `tests/cli/test_main.py` | Modify | Version assertion `0.5.1` |
| `docs/commands/tree-bi.md`, `.zh.md` | Modify | Stop signals, sidecar, SIGKILL, exit codes |
| `docs/commands/ai-integration.md`, `.zh.md` | Modify | `check_status` unknown + `stop_state` / `stop_state_error` |
| `skills/phyloai-workflow/SKILL.md` | Modify | Unknown handling, no false success |
| `skills/phyloai-workflow/references/error-catalog.md` | Modify | Exit 4 and 130/143/129 |
| `pyproject.toml`, `phyloai/__init__.py` | Modify | Version `0.5.1` (separate commit) |

---

## Task 0: Baseline and durable workspace

- [ ] Create worktree at `/home/chenwq/wsl/phyloai-pr1-fix` (not `/tmp`) from `36d4f96`.
- [ ] Run baseline; record base SHA, command, Python environment, result, and patch SHA-256.
- [ ] Keep a complete patch (tracked + untracked) in the review directory after each test-passing stage.

Report done — no commit.

---

## Task 1: `StopController` and signal lifecycle

**Files:** `phyloai/tree/bi.py`, `tests/tree/test_bi.py`.

- [ ] Add `StopController` with `first_signal`, `force_signal`, `stop_level`, `handled_level`, `phase`, latched commit signals.
- [ ] `requested`/`forced` derived; `pending_request = stop_level > handled_level`.
- [ ] `handled_level` means the stop action reached the **submission-ready condition**, not merely “request sent”.
- [ ] Install SIGINT/SIGTERM/SIGHUP; skip inherited `SIG_IGN`; restore only installed handlers.
- [ ] `install_handlers=False` still records signals (`note_signal` / `note_keyboard_interrupt`).
- [ ] Add/update Task 1 tests in `tests/tree/test_bi.py`; run them; do not defer to Task 6.

Report done — no commit.

---

## Task 2: Atomic JSON helper and `stop_state.json`

**Files:** `phyloai/core/schema.py`, `phyloai/tree/bi.py`, `tests/core/test_schema.py`, `tests/tree/test_bi.py`.

- [ ] Add `save_json_atomic(path, payload, *, fsync=False)`: same-directory temp file, flush/optional fsync, `os.replace`, cleanup on failure, preserve existing mode.
- [ ] `write_result_json()` uses the helper.
- [ ] `_write_stop_state()` writes minimal fields (`state`, `pid`, `requested_at`, `stop_signal`; forced adds `forced_at`/`force_signal`) with `fsync=True`.
- [ ] `_clear_stop_state()` helper; stale cleanup only after input/chain-command validation succeeds.
- [ ] Tests: serialization/replace failure preserves old target, temp cleanup, mode preservation, sidecar field preservation.

Report done — no commit.

---

## Task 3: Library shutdown, failed-chain cleanup, commit boundary

**Files:** `phyloai/tree/bi.py`, `tests/tree/test_bi.py`.

- [ ] Startup checks: check `stop.requested` before setup/launch and before each chain launch; do not launch remaining chains after a stop request.
- [ ] First request: minimal sidecar, soft stop, bounded `proc.wait(timeout=1)` loop, mark `handled_level=1` at submission-ready.
- [ ] Second request: forced sidecar, `terminate → bounded wait → kill → bounded wait`, skip final convergence, mark `handled_level=2`.
- [ ] Convergence tools use `Popen` + bounded polling; no blocking `subprocess.run()` when a controller is active.
- [ ] **Failed chain:** soft-stop remaining chains with the same bounded wait/controller checks; second signal forces them; preserve original chain-failure message.
- [ ] `begin_commit()`: phase first, then pending check; abort attempt and reset phase if pending; only `os.replace()` success is irreversible.
- [ ] Expose library finalization API for late requests; `run_bi_pb()` does not call `begin_commit()`.
- [ ] Direct library double Ctrl+C: first soft, second forced; do not put whole shutdown in one except body.
- [ ] Add/update Task 3 tests and run `tests/tree/test_bi.py`.

Report done — no commit.

---

## Task 4: CLI finalization, backup, business exit codes, sidecar cleanup

**Files:** `phyloai/cli/commands/tree.py`, `tests/cli/test_tree_bi.py`, `tests/cli/test_tree_bi_signals.py`.

- [ ] `bi_pb_command()` uses one `with StopController()` covering run, finalization, backup, publish, cleanup.
- [ ] Restore PR-before backup timing: backup only when a new result is ready to publish.
- [ ] Backup/serialize/flush/fsync/replace failure: preserve old result/sidecar; original success 0 → 4; forced 130/143/129 preserved; existing 1/2/3 unchanged.
- [ ] Final/error publication uses `fsync=True`; sidecar also `fsync=True`.
- [ ] **Publication exception isolation:** catch `(OSError, TypeError, ValueError)` from publication; circular-reference `ValueError` must not be classified as input error.
- [ ] **Business code first:** compute 0/2/forced before publication; publish failure must not turn tool failure 2 into 4.
- [ ] **Error result write failure:** `_write_error_result()` catches secondary write errors, warns, and preserves original 1/3.
- [ ] Late handoff: call library finalization API, update sidecar/payload, retry commit; never rerun analysis.
- [ ] Update CLI tests in this task; run `tests/cli/test_tree_bi.py tests/cli/test_tree_bi_signals.py`.

Report done — no commit.

---

## Task 5: `check_status` and Skill behavior

**Files:** `phyloai/mcp/tools/utils.py`, `tests/mcp/test_utils.py`, `skills/phyloai-workflow/SKILL.md`, `docs/commands/ai-integration.md`, `.zh.md`.

- [ ] Active PID → `running`; include `stop_state` or `stop_state_error`.
- [ ] Sidecar exists and not running → `unknown` + `stop_state`; include readable result but never treat as completion.
- [ ] Sidecar exists but unreadable/corrupt → `unknown` + `stop_state_error`; no stale-success fallback.
- [ ] No sidecar → existing result/job rules; no mtime/run identity.
- [ ] Tests: valid sidecar + old result, corrupt sidecar, running sidecar, sidecar removed.
- [ ] Skill/docs: unknown + stop_state/stop_state_error handling; no false success.

Report done — no commit.

---

## Task 6: Cross-phase tests and regression

- [ ] Second signal at stop entry, sidecar write, soft stop, wait, convergence, result assembly.
- [ ] Third/later signals do not repeat forced cleanup or change `force_signal`.
- [ ] Commit boundary: phase first, pending abort, commit-period signal does not change payload/exit code.
- [ ] Failed-chain branch: chain1 non-zero, chain2 ignores soft stop; second signal forces without blocking; original error preserved.
- [ ] Matrix-conversion failure preserves stale sidecar; chain-command validation failure also preserves it.
- [ ] Backup/serialize/flush/fsync/replace failure preserves old result/sidecar; original success 0 → 4; forced code preserved.
- [ ] Publication TypeError and circular-reference ValueError return 4; tool failure + publication failure returns 2; missing-tool error-write failure returns 3.
- [ ] Sidecar deletion failure warns and `check_status` remains unknown.
- [ ] Direct library double Ctrl+C; stop-before-launch does not start chains; internal `stop` absent from schema/params.
- [ ] Run six-file scoped suite and broader `tests/tree tests/cli tests/core tests/mcp`; record results.

Report done — no commit.

---

## Task 7: Spec and user-doc sync

- [ ] Verify tree-bi spec: stop signals, sidecar, result fields, exit codes, CLI/library boundary, failed-chain cleanup, stale-sidecar timing.
- [ ] Verify total design exit table has 4 and 128+signal.
- [ ] Verify JSON standard stop metadata and internal `stop`.
- [ ] Verify AI integration `check_status` sidecar/corrupt-sidecar rules.
- [ ] Update command docs and Skill/error catalog.

Report done — no commit.

---

## Task 8: Version 0.5.1 (separate commit)

- [ ] `pyproject.toml` and `phyloai/__init__.py` → `0.5.1`.
- [ ] `test_version_is_0_5_1`.
- [ ] Run `tests/cli/test_main.py -q`.

Report done — no commit.

---

## Acceptance Criteria

1. Normal run: success, `shutdown_state: completed`, no sidecar, exit 0.
2. Single SIGINT/SIGTERM/SIGHUP with chains exiting normally: soft stop, final convergence when eligible, success/0, `shutdown_state: stopped` + `stop_signal`, sidecar removed.
3. Second signal before commit: agreed direct-process termination/reap flow; convergence skipped; error payload `shutdown_state: forced` + `stop_signal`/`force_signal`; exit 130/143/129.
4. Hard kill after sidecar: `check_status` returns `unknown` + `stop_state`, never stale success.
5. SIGKILL before any stop record: documented pre-existing behavior; no startup marker, no run identity.
6. Commit failure: old result/sidecar preserved; original success 0 → 4; forced code preserved.
7. Library never writes `result.json`; CLI remains sole writer.
8. Atomic writes preserve old target and permissions.
9. Final/error result and sidecar use `fsync=True`.
10. Direct library Ctrl+C still soft-stops; second Ctrl+C forces.
11. All scoped tests pass; 6 readpb failures remain separate.
