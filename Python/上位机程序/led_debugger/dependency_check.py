from __future__ import annotations

import re
import locale
import subprocess
import sys
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path


@dataclass(frozen=True)
class RequirementItem:
    """requirements.txt 中的一条依赖要求。"""

    name: str
    operator: str | None
    required_version: str | None
    raw_line: str


@dataclass(frozen=True)
class MissingRequirement:
    """缺失或版本不满足的依赖。"""

    requirement: RequirementItem
    reason: str


def parse_requirements(requirements_path: Path) -> list[RequirementItem]:
    """解析 requirements.txt 中当前项目使用到的简单依赖格式。"""
    items: list[RequirementItem] = []
    if not requirements_path.exists():
        return items

    for raw_line in requirements_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        match = re.match(r"^([A-Za-z0-9_.-]+)\s*(>=|==)?\s*([A-Za-z0-9_.+-]+)?", line)
        if not match:
            continue
        name, operator, required_version = match.groups()
        items.append(
            RequirementItem(
                name=name,
                operator=operator,
                required_version=required_version,
                raw_line=line,
            )
        )
    return items


def find_missing_requirements(requirements_path: Path) -> list[MissingRequirement]:
    """检查 requirements.txt 中依赖是否已安装并满足最低版本。"""
    missing: list[MissingRequirement] = []
    for item in parse_requirements(requirements_path):
        try:
            installed_version = version(item.name)
        except PackageNotFoundError:
            missing.append(MissingRequirement(item, "未安装"))
            continue

        if item.operator == ">=" and item.required_version is not None:
            if _version_tuple(installed_version) < _version_tuple(item.required_version):
                missing.append(
                    MissingRequirement(
                        item,
                        f"版本过低，当前 {installed_version}，需要 >= {item.required_version}",
                    )
                )
        elif item.operator == "==" and item.required_version is not None:
            if installed_version != item.required_version:
                missing.append(
                    MissingRequirement(
                        item,
                        f"版本不一致，当前 {installed_version}，需要 == {item.required_version}",
                    )
                )
    return missing


def install_requirements(requirements_path: Path) -> tuple[bool, str]:
    """安装 requirements.txt，返回是否成功和 pip 输出。"""
    command = [sys.executable, "-m", "pip", "install", "-r", str(requirements_path)]
    result = subprocess.run(command, capture_output=True, text=True)
    output = "\n".join(part for part in (result.stdout, result.stderr) if part)
    return result.returncode == 0, output.strip()


def decode_process_output(data: bytes) -> str:
    """解码 pip 子进程输出，兼容 Windows 中文路径和不同控制台编码。"""
    candidate_encodings = [
        "utf-8",
        locale.getpreferredencoding(False),
        "gbk",
    ]
    tried: set[str] = set()
    for encoding in candidate_encodings:
        normalized = encoding.lower()
        if normalized in tried:
            continue
        tried.add(normalized)
        try:
            return data.decode(encoding)
        except UnicodeDecodeError:
            continue
    return data.decode("utf-8", errors="replace")


def show_dependency_dialog(missing: list[MissingRequirement], requirements_path: Path) -> bool:
    """显示缺少依赖弹窗；用户要求继续启动返回 True，退出返回 False。

    PySide6 缺失时无法创建该弹窗，调用方需要先判断并退化到命令行提示。
    """
    from PySide6.QtCore import QProcess, QProcessEnvironment
    from PySide6.QtGui import QColor, QFont, QLinearGradient, QPainter, QPainterPath, QPen, QTextCursor
    from PySide6.QtWidgets import (
        QApplication,
        QDialog,
        QHBoxLayout,
        QLabel,
        QPushButton,
        QTextEdit,
        QVBoxLayout,
        QWidget,
    )

    class DependencyDialog(QDialog):
        def __init__(self) -> None:
            super().__init__()
            self.should_continue = False
            self.install_process: QProcess | None = None
            self.setWindowTitle("依赖检查")
            self.setModal(True)
            self.resize(520, 360)
            self._build_ui()

        def _build_ui(self) -> None:
            layout = QVBoxLayout(self)
            layout.setContentsMargins(22, 22, 22, 22)
            layout.setSpacing(14)

            title = QLabel("程序缺少必要依赖")
            title.setObjectName("dialogTitle")

            message = QLabel("检测到当前 Python 环境未满足 requirements.txt，请安装后继续启动。")
            message.setWordWrap(True)
            message.setObjectName("dialogMessage")

            self.detail = QTextEdit()
            self.detail.setObjectName("detailBox")
            self.detail.setReadOnly(True)
            self.detail.setPlainText(
                "缺少依赖：\n"
                + "\n".join(f"- {item.requirement.raw_line}：{item.reason}" for item in missing)
            )

            self.status = QLabel("")
            self.status.setObjectName("dialogStatus")
            self.status.setWordWrap(True)

            buttons = QHBoxLayout()
            buttons.addStretch(1)
            self.install_button = QPushButton("自动安装")
            self.install_button.setObjectName("primaryButton")
            self.exit_button = QPushButton("退出程序")
            self.exit_button.setObjectName("secondaryButton")
            self.install_button.clicked.connect(self._start_install)
            self.exit_button.clicked.connect(self.reject)
            buttons.addWidget(self.install_button)
            buttons.addWidget(self.exit_button)

            layout.addWidget(title)
            layout.addWidget(message)
            layout.addWidget(self.detail, 1)
            layout.addWidget(self.status)
            layout.addLayout(buttons)
            self.setStyleSheet(
                """
                QDialog {
                    background: qlineargradient(
                        x1: 0, y1: 0, x2: 1, y2: 1,
                        stop: 0 #f8fafc,
                        stop: 1 #e8eef7
                    );
                    color: #0f172a;
                    font-family: "Microsoft YaHei UI", "Segoe UI", sans-serif;
                    font-size: 14px;
                }
                #dialogTitle {
                    font-size: 20px;
                    font-weight: 600;
                    color: #0f172a;
                }
                #dialogMessage, #dialogStatus {
                    color: #475569;
                }
                #detailBox {
                    background: rgba(255, 255, 255, 0.62);
                    border: 1px solid rgba(255, 255, 255, 0.78);
                    border-radius: 14px;
                    padding: 10px;
                    color: #0f172a;
                }
                QPushButton {
                    min-width: 92px;
                    min-height: 34px;
                    border-radius: 12px;
                    padding: 4px 14px;
                }
                #primaryButton {
                    background: rgba(219, 234, 254, 0.76);
                    border: 1px solid rgba(59, 130, 246, 0.72);
                    color: #0f172a;
                }
                #secondaryButton {
                    background: rgba(255, 255, 255, 0.62);
                    border: 1px solid rgba(203, 213, 225, 0.82);
                    color: #0f172a;
                }
                QPushButton:hover {
                    background: rgba(255, 255, 255, 0.86);
                }
                """
            )

        def _start_install(self) -> None:
            self.install_button.setEnabled(False)
            self.exit_button.setEnabled(False)
            self.status.setText("正在安装依赖，请稍候...")
            self.detail.clear()

            self.install_process = QProcess(self)
            self.install_process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
            # 强制 pip 子进程用 UTF-8 输出，避免中文项目路径显示成乱码。
            environment = QProcessEnvironment.systemEnvironment()
            environment.insert("PYTHONUTF8", "1")
            environment.insert("PYTHONIOENCODING", "utf-8")
            self.install_process.setProcessEnvironment(environment)
            self.install_process.readyReadStandardOutput.connect(self._append_install_output)
            self.install_process.finished.connect(self._finish_install)
            self.install_process.start(
                sys.executable,
                ["-m", "pip", "install", "-r", str(requirements_path)],
            )

        def _append_install_output(self) -> None:
            if self.install_process is None:
                return
            output = decode_process_output(bytes(self.install_process.readAllStandardOutput()))
            if not output:
                return
            self.detail.moveCursor(QTextCursor.MoveOperation.End)
            self.detail.insertPlainText(output)
            self.detail.moveCursor(QTextCursor.MoveOperation.End)

        def _finish_install(self, exit_code: int, _exit_status) -> None:
            self._append_install_output()
            self.exit_button.setEnabled(True)
            if exit_code == 0:
                self.status.setText("安装完成，请重新启动程序。")
                self.install_button.setText("安装完成")
                self.install_button.setEnabled(False)
                self.exit_button.setText("退出程序")
            else:
                self.status.setText("依赖安装失败，请检查网络或 pip 环境。")
                self.install_button.setEnabled(True)
                self.exit_button.setEnabled(True)

    app = QApplication.instance() or QApplication(sys.argv)
    dialog = DependencyDialog()
    result = dialog.exec()
    return result == QDialog.DialogCode.Accepted and dialog.should_continue


def _version_tuple(value: str) -> tuple[int, ...]:
    """把常见版本号转成可比较元组，满足当前 requirements 的 >= 检查。"""
    parts = re.findall(r"\d+", value)
    return tuple(int(part) for part in parts)
