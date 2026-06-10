from __future__ import annotations

from pathlib import Path

from scripts import entry_validate
from scripts import run_validation_suite


def test_entry_validate_invokes_validation_script(monkeypatch, capsys):
    calls: list[list[str]] = []

    class DummyProc:
        returncode = 0

    def fake_run(cmd):
        calls.append(cmd)
        return DummyProc()

    monkeypatch.setattr(entry_validate.subprocess, "run", fake_run)

    rc = entry_validate.main()

    assert rc == 0
    assert calls == [[entry_validate.sys.executable, "-m", "scripts.run_validation_suite"]]
    out = capsys.readouterr().out
    assert "[entry_validate] Starting validation" in out
    assert "[entry_validate] Validation finished with rc=0" in out


def test_run_validation_suite_replaces_stale_report(monkeypatch, tmp_path, capsys):
    stale_report = tmp_path / "validation_suite.md"
    stale_report.write_text("STALE PASS", encoding="utf-8")

    def fake_run(cmd):
        return 0, f"completed {' '.join(cmd[-2:])}"

    monkeypatch.setattr(run_validation_suite, "_run", fake_run)

    rc = run_validation_suite.main(["--output_md", str(stale_report)])

    assert rc == 0
    text = stale_report.read_text(encoding="utf-8")
    assert "STALE PASS" not in text
    assert "Overall command status: PASS" in text
    assert "falsifiability_check [empirical]: PASS" in text
    assert not Path(str(stale_report) + ".tmp").exists()
    out = capsys.readouterr().out
    assert "[validation] Step 1/" in out
    assert "PASS falsifiability_check (rc=0)" in out
