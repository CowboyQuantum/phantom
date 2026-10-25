import phantom


def test_imports_and_version():
    assert phantom.__version__
    assert isinstance(phantom.MODEL_VERSION, int)
