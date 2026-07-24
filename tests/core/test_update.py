"""Tests for phyloai.core.update."""

from __future__ import annotations

import json
from unittest.mock import patch

import pytest


class TestSemVer:
    def test_parse_simple(self) -> None:
        from phyloai.core.update import _SemVer

        v = _SemVer.parse("1.2.3")
        assert v.major == 1
        assert v.minor == 2
        assert v.patch == 3
        assert v.prerelease == ()
        assert v.build == ""

    def test_parse_with_v_prefix(self) -> None:
        from phyloai.core.update import _SemVer

        v = _SemVer.parse("v0.3.0")
        assert v.major == 0
        assert v.minor == 3
        assert v.patch == 0

    def test_parse_with_prerelease(self) -> None:
        from phyloai.core.update import _SemVer

        v = _SemVer.parse("1.0.0-alpha.1")
        assert v.prerelease == ("alpha", 1)

    def test_parse_with_build(self) -> None:
        from phyloai.core.update import _SemVer

        v = _SemVer.parse("1.0.0+build123")
        assert v.build == "build123"
        assert v.prerelease == ()

    def test_ordering_simple(self) -> None:
        from phyloai.core.update import _SemVer

        assert _SemVer.parse("1.0.0") < _SemVer.parse("2.0.0")
        assert _SemVer.parse("1.0.0") < _SemVer.parse("1.1.0")
        assert _SemVer.parse("1.0.0") < _SemVer.parse("1.0.1")

    def test_ordering_prerelease(self) -> None:
        from phyloai.core.update import _SemVer

        assert _SemVer.parse("1.0.0-alpha") < _SemVer.parse("1.0.0")
        assert _SemVer.parse("1.0.0-alpha") < _SemVer.parse("1.0.0-beta")
        assert _SemVer.parse("1.0.0-alpha.1") < _SemVer.parse("1.0.0-alpha.2")
        assert _SemVer.parse("1.0.0-alpha.2") < _SemVer.parse("1.0.0-alpha.10")
        assert _SemVer.parse("1.0.0-rc.2") < _SemVer.parse("1.0.0-rc.10")

    def test_build_does_not_affect_order(self) -> None:
        from phyloai.core.update import _SemVer

        a = _SemVer.parse("1.0.0+build1")
        b = _SemVer.parse("1.0.0+build2")
        assert not (a < b)
        assert not (b < a)

    def test_v_prefix_tags_compare(self) -> None:
        from phyloai.core.update import _SemVer

        assert _SemVer.parse("v0.2.0") < _SemVer.parse("v0.3.0")

    def test_not_equal(self) -> None:
        from phyloai.core.update import _SemVer

        assert not (_SemVer.parse("0.3.0") < _SemVer.parse("0.3.0"))


class TestFetchLatest:
    def test_returns_none_on_network_error(self) -> None:
        from phyloai.core.update import _fetch_latest_version
        from urllib.error import URLError

        with patch("phyloai.core.update.urlopen", side_effect=URLError("fail")):
            assert _fetch_latest_version() is None

    def test_returns_highest_semver_tag(self) -> None:
        from phyloai.core.update import _fetch_latest_version

        mock_response = json.dumps([
            {"name": "v0.1.0"},
            {"name": "v0.3.0"},
            {"name": "v0.2.0"},
            {"name": "nightly"},
        ]).encode()

        class FakeResponse:
            headers = {"Link": ""}

            @staticmethod
            def read():
                return mock_response

            def __enter__(self):
                return self

            def __exit__(self, *args):
                pass

        with patch("phyloai.core.update.urlopen", return_value=FakeResponse()):
            result = _fetch_latest_version()
        assert result == "v0.3.0"

    def test_returns_none_for_empty_tags(self) -> None:
        from phyloai.core.update import _fetch_latest_version

        class FakeResponse:
            headers = {"Link": ""}

            @staticmethod
            def read():
                return b"[]"

            def __enter__(self):
                return self

            def __exit__(self, *args):
                pass

        with patch("phyloai.core.update.urlopen", return_value=FakeResponse()):
            assert _fetch_latest_version() is None

    def test_pagination_follows_multiple_pages(self) -> None:
        from phyloai.core.update import _fetch_latest_version

        page1 = json.dumps([{"name": "v0.1.0"}, {"name": "v0.2.0"}]).encode()
        page2 = json.dumps([{"name": "v0.3.0"}, {"name": "v0.4.0"}]).encode()

        call_count = [0]

        class FakeResponse:
            def __init__(self, data, link=""):
                self._data = data
                self.headers = {"Link": link}

            def read(self):
                return self._data

            def __enter__(self):
                return self

            def __exit__(self, *args):
                pass

        def fake_urlopen(req, timeout=None):
            call_count[0] += 1
            url = req.full_url if hasattr(req, "full_url") else str(req)
            if call_count[0] == 1:
                return FakeResponse(page1, link='<https://api.github.com/repos/xtmtd/phyloai/tags?page=2>; rel="next"')
            else:
                return FakeResponse(page2)

        with patch("phyloai.core.update.urlopen", side_effect=fake_urlopen):
            result = _fetch_latest_version()
        assert result == "v0.4.0"
        assert call_count[0] == 2


class TestCheckUpdate:
    def test_up_to_date(self) -> None:
        from phyloai.core.update import check_update

        with patch("phyloai.core.update._fetch_latest_version", return_value="v0.3.0"):
            with patch("phyloai.core.update.__version__", "0.3.0"):
                result = check_update()
        assert result["status"] == "up_to_date"

    def test_update_available(self) -> None:
        from phyloai.core.update import check_update

        with patch("phyloai.core.update._fetch_latest_version", return_value="v0.4.0"):
            with patch("phyloai.core.update.__version__", "0.3.0"):
                result = check_update()
        assert result["status"] == "available"
        assert result["current"] == "0.3.0"
        assert result["latest"] == "v0.4.0"

    def test_fetch_failure(self) -> None:
        from phyloai.core.update import check_update

        with patch("phyloai.core.update._fetch_latest_version", return_value=None):
            result = check_update()
        assert result["status"] == "error"


class TestRunUpdate:
    def test_up_to_date_exits_0(self) -> None:
        from phyloai.core.update import run_update

        with patch("phyloai.core.update.check_update", return_value={"status": "up_to_date", "current": "0.3.0", "latest": "v0.3.0"}):
            assert run_update(confirm=False) == 0

    def test_fetch_error_exits_1(self) -> None:
        from phyloai.core.update import run_update

        with patch("phyloai.core.update.check_update", return_value={"status": "error", "message": "fail"}):
            assert run_update(confirm=False) == 1

    def test_confirm_yes_installs(self) -> None:
        from phyloai.core.update import run_update

        with patch("phyloai.core.update.check_update", return_value={"status": "available", "current": "0.3.0", "latest": "v0.4.0"}):
            with patch("phyloai.core.update.subprocess.run", return_value=type("Fake", (), {"returncode": 0})()) as mock_run:
                result = run_update(confirm=True)
        assert result == 0
        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert "@v0.4.0" in cmd[-1]
