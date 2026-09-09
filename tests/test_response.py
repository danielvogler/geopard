"""The object every match returns.

``dtw_match`` never raises for a match that simply did not happen -- it returns
a response carrying a flag and, when relevant, an error string. Callers branch
on that, so its shape is part of the public contract.
"""

from __future__ import annotations

from datetime import timedelta

import pytest

from geopard import GeopardResponse


def _response(**overrides):
    fields = {
        "time": timedelta(seconds=600),
        "dtw": 0.1,
        "start_point": [47.1, 9.1, 1000.0, None],
        "end_point": [47.2, 9.2, 1100.0, None],
        "match_flag": 2,
        "error": None,
    }
    fields.update(overrides)
    return GeopardResponse(**fields)


def test_carries_every_documented_attribute():
    response = _response()

    assert response.time == timedelta(seconds=600)
    assert response.dtw == 0.1
    assert response.start_point[0] == 47.1
    assert response.end_point[0] == 47.2
    assert response.match_flag == 2
    assert response.error is None


@pytest.mark.parametrize("flag", [1, 2])
def test_positive_flags_are_successes(flag):
    assert _response(match_flag=flag).is_success() is True


@pytest.mark.parametrize("flag", [0, -1, -2])
def test_zero_and_negative_flags_are_not(flag):
    assert _response(match_flag=flag).is_success() is False


def test_error_defaults_to_none():
    response = GeopardResponse(None, None, None, None, -1)

    assert response.error is None


def test_accepts_positional_arguments_in_the_documented_order():
    """The documented positional order must not drift.

    The matcher constructs responses positionally, and so does downstream code.
    """
    response = GeopardResponse(timedelta(seconds=1), 0.5, [1], [2], -2, "boom")

    assert response.time == timedelta(seconds=1)
    assert response.dtw == 0.5
    assert response.start_point == [1]
    assert response.end_point == [2]
    assert response.match_flag == -2
    assert response.error == "boom"


def test_is_success_is_callable_not_a_bare_bool():
    """The shipped examples call it. Turning it into a property would break them."""
    assert callable(_response().is_success)


def test_parse_response_logs_without_raising(gp, caplog):
    """A convenience printer; it must survive a fully populated response."""
    import logging

    logging.disable(logging.NOTSET)
    with caplog.at_level(logging.INFO):
        gp.parse_response(_response())

    assert "Response DTW" in caplog.text
