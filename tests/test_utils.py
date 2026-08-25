import sys
from pathlib import Path

import pandas as pd
import pytest

from qpcr.utils import gradient_styles, rdylgn_hex

REPO = Path(__file__).resolve().parents[1]


class TestRdYlGnGradient:
    """The results-table gradient, which must not need matplotlib.

    `Styler.background_gradient` imports matplotlib at CALL time. The pin was
    dropped as unused because nothing writes `import matplotlib`, and every
    styled results table then raised "background_gradient requires matplotlib"
    — but only where matplotlib was genuinely absent, so a dev box with a
    leftover copy stayed green while CI and Streamlit Cloud broke.
    """

    def test_ramp_runs_red_through_yellow_to_green(self):
        assert rdylgn_hex(0.0) == "#a50026"      # dark red
        assert rdylgn_hex(0.5) == "#ffffbf"      # pale yellow
        assert rdylgn_hex(1.0) == "#006837"      # dark green

    def test_ramp_clamps_outside_the_unit_interval(self):
        assert rdylgn_hex(-5.0) == rdylgn_hex(0.0)
        assert rdylgn_hex(99.0) == rdylgn_hex(1.0)

    def test_values_map_across_the_configured_range(self):
        # vmin=0, vmax=3 is what the results table uses: 1.5 is mid-ramp.
        styles = gradient_styles([0.0, 1.5, 3.0], vmin=0, vmax=3)
        assert "#a50026" in styles[0]
        assert "#ffffbf" in styles[1]
        assert "#006837" in styles[2]

    def test_out_of_range_values_clamp_rather_than_wrap(self):
        styles = gradient_styles([-2.0, 12.0], vmin=0, vmax=3)
        assert "#a50026" in styles[0]
        assert "#006837" in styles[1]

    def test_missing_values_are_left_unpainted(self):
        """An unmeasurable cell must not be coloured as if it sat on the scale."""
        styles = gradient_styles([float("nan"), None, "n/a"], vmin=0, vmax=3)
        assert styles == ["", "", ""]

    def test_dark_cells_get_light_text(self):
        dark = gradient_styles([0.0], vmin=0, vmax=3)[0]
        mid = gradient_styles([1.5], vmin=0, vmax=3)[0]
        assert "color: #f1f1f1" in dark      # dark red background
        assert "color: #000000" in mid       # pale yellow background

    def test_zero_width_range_does_not_divide_by_zero(self):
        assert gradient_styles([1.0], vmin=2.0, vmax=2.0) != [""]

    def test_styles_the_results_table_with_matplotlib_absent(self, monkeypatch):
        """The actual regression: reproduce a clean install and render.

        Setting the module to None makes `import matplotlib` raise ImportError,
        which is precisely what CI hit.
        """
        monkeypatch.setitem(sys.modules, "matplotlib", None)
        with pytest.raises(ImportError):
            import matplotlib  # noqa: F401

        df = pd.DataFrame({"Condition": ["Ctrl", "Test"],
                           "Fold_Change": [1.0, 2.4]})
        styled = df.style.apply(gradient_styles, subset=["Fold_Change"],
                                vmin=0, vmax=3)
        styled._compute()   # what st.dataframe triggers, and where it blew up
        assert styled.to_html()

    def test_app_does_not_reintroduce_background_gradient(self):
        """Guard the whole class of bug, not just the one call site.

        Comment lines are skipped so the notes explaining *why* this is banned
        do not themselves trip the guard.
        """
        offenders = []
        for path in [REPO / "streamlit qpcr analysis v1.py",
                     *(REPO / "qpcr").rglob("*.py")]:
            for lineno, line in enumerate(
                path.read_text(encoding="utf-8").splitlines(), 1
            ):
                if line.lstrip().startswith("#"):
                    continue
                if ".background_gradient(" in line:
                    offenders.append(f"{path.name}:{lineno}")
        assert not offenders, (
            f"{offenders} call Styler.background_gradient, which imports "
            "matplotlib at call time; matplotlib is not in requirements.txt"
        )

    def test_matplotlib_is_not_a_dependency(self):
        """If this pin ever comes back, the workaround above is redundant."""
        pins = (REPO / "requirements.txt").read_text(encoding="utf-8").lower()
        assert "matplotlib" not in pins
        assert "seaborn" not in pins


class TestNaturalSortKey:
    def test_natural_sort_key_basic_numbers(self, mock_streamlit):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        natural_sort_key = spec.natural_sort_key
        
        samples = ['Sample10', 'Sample2', 'Sample1', 'Sample20']
        sorted_samples = sorted(samples, key=natural_sort_key)
        
        assert sorted_samples == ['Sample1', 'Sample2', 'Sample10', 'Sample20']

    def test_natural_sort_key_mixed_content(self, mock_streamlit):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        natural_sort_key = spec.natural_sort_key
        
        samples = ['A10B2', 'A2B1', 'A1B10']
        sorted_samples = sorted(samples, key=natural_sort_key)
        
        assert sorted_samples == ['A1B10', 'A2B1', 'A10B2']

    def test_natural_sort_key_no_numbers(self, mock_streamlit):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        natural_sort_key = spec.natural_sort_key
        
        samples = ['Charlie', 'Alpha', 'Bravo']
        sorted_samples = sorted(samples, key=natural_sort_key)
        
        assert sorted_samples == ['Alpha', 'Bravo', 'Charlie']

    def test_natural_sort_key_case_insensitive(self, mock_streamlit):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        natural_sort_key = spec.natural_sort_key
        
        samples = ['sample1', 'SAMPLE2', 'Sample3']
        sorted_samples = sorted(samples, key=natural_sort_key)
        
        assert sorted_samples == ['sample1', 'SAMPLE2', 'Sample3']

    def test_natural_sort_key_empty_string(self, mock_streamlit):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        natural_sort_key = spec.natural_sort_key
        
        result = natural_sort_key('')
        assert result == ['']

    def test_natural_sort_key_none_handling(self, mock_streamlit):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        natural_sort_key = spec.natural_sort_key
        
        result = natural_sort_key(None)
        assert isinstance(result, list)
