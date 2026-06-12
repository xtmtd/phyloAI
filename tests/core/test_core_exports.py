def test_core_exposes_checkpoint_helpers() -> None:
    from phyloai.core import (
        Checkpoint,
        CheckpointTask,
        canonical_params_hash,
        load_checkpoint,
        save_checkpoint_atomic,
        summarize_resume_tasks,
        validate_resume_params,
    )

    assert Checkpoint is not None
    assert CheckpointTask is not None
    assert callable(canonical_params_hash)
    assert callable(load_checkpoint)
    assert callable(save_checkpoint_atomic)
    assert callable(summarize_resume_tasks)
    assert callable(validate_resume_params)
