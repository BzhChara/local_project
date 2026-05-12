from pathlib import Path
import sys

from led_debugger.dependency_check import (
    MissingRequirement,
    RequirementItem,
    find_missing_requirements,
    show_dependency_dialog,
)


def _show_dependency_check_demo(requirements_path: Path) -> int:
    """人工预览依赖检查弹窗，不参与正常启动流程。"""
    demo_missing = [
        MissingRequirement(
            RequirementItem("pyserial", ">=", "3.5", "pyserial>=3.5"),
            "演示：未安装",
        ),
        MissingRequirement(
            RequirementItem("pyqtgraph", ">=", "0.13", "pyqtgraph>=0.13"),
            "演示：版本过低，当前 0.12，需要 >= 0.13",  
        ),
    ]
    return 0 if show_dependency_dialog(demo_missing, requirements_path) else 1


if __name__ == "__main__":
    # 先检查 requirements.txt，避免主窗口导入阶段因缺依赖直接崩溃。
    project_root = Path(__file__).resolve().parent
    requirements_path = project_root / "requirements.txt"
    if "--show-dependency-check-demo" in sys.argv:
        sys.exit(_show_dependency_check_demo(requirements_path))

    missing = find_missing_requirements(requirements_path)
    if missing:
        py_side_missing = any(item.requirement.name.lower() == "pyside6" for item in missing)
        if py_side_missing:
            print("缺少 PySide6，无法显示图形化依赖检查窗口。")
            print(f"请执行：{sys.executable} -m pip install -r {requirements_path}")
            sys.exit(1)

        if not show_dependency_dialog(missing, requirements_path):
            sys.exit(1)

        missing = find_missing_requirements(requirements_path)
        if missing:
            print("依赖仍未满足：")
            for item in missing:
                print(f"- {item.requirement.raw_line}：{item.reason}")
            sys.exit(1)

    from led_debugger.main_window import run

    sys.exit(run())
