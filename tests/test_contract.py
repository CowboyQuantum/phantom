"""P0 for Step 8: the consumer contract (`python -m phantom.contract`) holds on the reference device and on
other draws; the module's exit status reports it."""

import pytest

from phantom.contract import CHECKS, main, reference_device, run


@pytest.mark.parametrize("seed", [2, 4, 7])
def test_every_contract_check_holds(seed: int):
    results = run(reference_device(seed))
    assert len(results) == len(CHECKS)
    failed = [f"{c.name}: {c.detail}" for c in results if not c.ok]
    assert not failed, failed
    assert not [c.name for c in results if c.detail.startswith("skipped")], "every check ran"


def test_the_module_exits_zero_and_reads_device_files(tmp_path, capsys):
    path = tmp_path / "device.json"
    reference_device(5).to_json(path)
    assert main([str(path)]) == 0
    out = capsys.readouterr().out
    assert "contract: OK" in out and str(path) in out


def test_a_broken_promise_is_reported_not_raised(monkeypatch):
    from phantom.device import Device

    monkeypatch.setattr(Device, "knowable", lambda self: None)
    results = {c.name: c for c in run(reference_device())}
    assert not results["knowable"].ok and "AttributeError" in results["knowable"].detail
    assert results["truth"].ok, "the other checks still ran"
