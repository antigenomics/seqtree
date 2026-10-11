"""Setup combines extras and never reports success after a failed install."""
import os
from pathlib import Path
import subprocess

import pytest


@pytest.mark.skipif(os.name == "nt", reason="POSIX setup script")
def test_setup_extras_and_install_failure(tmp_path):
    uv = tmp_path / "uv"
    uv.write_text('#!/bin/sh\nprintf "%s\\n" "$*"\nexit "${TEST_UV_STATUS:-0}"\n')
    uv.chmod(0o755)
    env = dict(os.environ, PATH=str(tmp_path) + os.pathsep + os.environ["PATH"])
    root = Path(__file__).resolve().parents[2]
    ok = subprocess.run(["sh", "setup.sh", "--tests", "--bench"], cwd=root,
                        env=env, capture_output=True, text=True, timeout=10)
    assert ok.returncode == 0 and "pip install -e .[test,bench]" in ok.stdout
    failed = subprocess.run(["sh", "setup.sh", "--tests"], cwd=root,
                            env=dict(env, TEST_UV_STATUS="42"), capture_output=True,
                            text=True, timeout=10)
    assert failed.returncode == 42 and "Done." not in failed.stdout
