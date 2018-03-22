import pytest

def pytest_addoption(parser):
    parser.addoption("--comparison_dir", action="store", default=None, help="The path to the comparison directory")

@pytest.fixture(name='comparison_dir')
def fixture_comparison_dir(request):
    return request.config.option.comparison_dir