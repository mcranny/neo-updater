"""Zero-activation entrypoint for Asteroid Intercept Planner."""

from __future__ import annotations

import importlib.util
import os
import subprocess
import sys
import venv
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RUNTIME = ROOT / ".venv"


def _runtime_python() -> Path:
    if os.name == "nt":
        return RUNTIME / "Scripts" / "python.exe"
    return RUNTIME / "bin" / "python"


def _desktop_dependencies_available() -> bool:
    modules = (
        "PySide6",
        "pyqtgraph",
        "OpenGL",
        "flask",
        "numpy",
    )
    return all(importlib.util.find_spec(module) is not None for module in modules)


def _install_runtime(python: Path) -> None:
    print("Preparing the Asteroid Intercept Planner runtime…", flush=True)
    subprocess.check_call([str(python), "-m", "pip", "install", "--upgrade", "pip"])
    subprocess.check_call([str(python), "-m", "pip", "install", "-e", f"{ROOT}[desktop]"])


def _ensure_runtime() -> None:
    if sys.version_info < (3, 11):  # noqa: UP036 - bootstrap must fail clearly
        raise SystemExit("Python 3.11 or newer is required.")
    runtime_python = _runtime_python()
    if Path(sys.executable).resolve() != runtime_python.resolve():
        if not runtime_python.exists():
            print("First launch: creating a private local Python environment…", flush=True)
            venv.EnvBuilder(with_pip=True).create(RUNTIME)
        os.execv(
            str(runtime_python),
            [str(runtime_python), str(Path(__file__).resolve()), *sys.argv[1:]],
        )
    if not _desktop_dependencies_available():
        _install_runtime(runtime_python)


def main() -> int:
    _ensure_runtime()
    from app.integrated_app import main as launch_application

    return launch_application()


if __name__ == "__main__":
    raise SystemExit(main())
