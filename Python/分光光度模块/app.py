import sys


def main():
    try:
        from spectro_ui.main_window import run
    except ModuleNotFoundError as exc:
        if exc.name == "PySide6":
            print("缺少 PySide6，无法启动图形界面。")
            print("请执行：python -m pip install -r requirements.txt")
            return 1
        raise

    return run()


if __name__ == "__main__":
    sys.exit(main())
