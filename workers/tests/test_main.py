import inspect
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import main


def test_worker_entrypoint_is_async():
    assert inspect.iscoroutinefunction(main.main)
