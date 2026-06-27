# Error Catalog

Use this for exit code 1 and 3 errors. Exit code 2 is external-tool stderr diagnosis.

## Exit 1: Missing Required Parameter
Pattern: `Missing option`, `is required`
Fix: Re-render the parameter card and request the missing field.

## Exit 1: Output Directory Exists
Pattern: `Output directory exists and is not empty`
Fix: Use a new `--output-dir`, or add `--overwrite` after user confirmation.

## Exit 1: Invalid Enum
Pattern: `Invalid value for`, `Choose from`
Fix: Show allowed enum values from runtime schema.

## Exit 1: File Not Found
Pattern: `does not exist`, `not found`
Fix: Ask user to confirm the path; prefer absolute paths for recovery.

## Exit 1: Empty File
Pattern: `empty`
Fix: Exclude the file or regenerate upstream input.

## Exit 1: Unsupported Format
Pattern: `Unsupported`, `unrecognized format`
Fix: Run `pretree convert` or set `--input-format` if appropriate.

## Exit 1: Overwrite Resume Conflict
Pattern: `--overwrite and --resume are mutually exclusive`
Fix: Choose one: overwrite for fresh run, resume for interrupted run.

## Exit 1: Blocked Tool Args
Pattern: `blocked flag`, `managed by PhyloAI`
Fix: Remove path/control flags from `--tool-args`; use PhyloAI parameters instead.

## Exit 3: Required Tool Missing
Pattern: `not installed`, `not detectable`, `not found`
Fix: Run `doctor`; install the missing tool or configure its path.

## Exit 3: Runtime Missing
Pattern: `Java`, `Julia`, `MPI`
Fix: Install the required runtime and rerun `doctor`.
