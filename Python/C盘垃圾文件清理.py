from __future__ import annotations

import argparse
import logging
import os
import shutil
import stat
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable


PROTECTED_PREFIXES = (
    Path(r"C:\Windows\System32"),
    Path(r"C:\Windows\SysWOW64"),
    Path(r"C:\Program Files"),
    Path(r"C:\Program Files (x86)"),
)

GLOB_PATTERNS = ("*.tmp", "*.log")


@dataclass(slots=True)
class ScanTarget:
    name: str
    path: Path
    mode: str  # "tree" or "glob"


@dataclass(slots=True)
class FileCandidate:
    path: Path
    size: int
    source: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="扫描并清理 Windows 常见垃圾文件。默认 dry-run，仅输出结果。"
    )
    parser.add_argument(
        "--clean",
        "--force",
        dest="clean",
        action="store_true",
        help="实际执行删除。默认只做 dry-run。",
    )
    parser.add_argument(
        "--exclude",
        action="append",
        default=[],
        metavar="PATH",
        help="排除目录或文件，可重复传入。",
    )
    parser.add_argument(
        "--log-file",
        default=str(default_log_file()),
        help="日志文件路径，默认写到脚本同目录下 logs。",
    )
    return parser.parse_args()


def default_log_file() -> Path:
    base_dir = Path(__file__).resolve().parent / "logs"
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return base_dir / f"cleanup_{timestamp}.log"


def setup_logging(log_file: str) -> Path:
    log_path = Path(log_file).expanduser().resolve()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_path, encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )
    return log_path


def build_targets() -> list[ScanTarget]:
    local_app_data = Path(os.environ.get("LOCALAPPDATA", r"C:\Users\Default\AppData\Local"))
    user_profile = Path(os.environ.get("USERPROFILE", r"C:\Users\Default"))
    windir = Path(os.environ.get("WINDIR", r"C:\Windows"))
    system_drive = Path(os.environ.get("SystemDrive", "C:") + "\\")

    return [
        ScanTarget("Temp 临时文件", local_app_data / "Temp", "tree"),
        ScanTarget("Windows Temp", windir / "Temp", "tree"),
        ScanTarget("Windows 更新缓存", windir / "SoftwareDistribution" / "Download", "tree"),
        ScanTarget("Chrome 缓存", local_app_data / "Google" / "Chrome" / "User Data" / "Default" / "Cache", "tree"),
        ScanTarget("Edge 缓存", local_app_data / "Microsoft" / "Edge" / "User Data" / "Default" / "Cache", "tree"),
        ScanTarget("Firefox 缓存", local_app_data / "Mozilla" / "Firefox" / "Profiles", "tree"),
        ScanTarget("回收站", system_drive / "$Recycle.Bin", "tree"),
        ScanTarget("用户下载目录中的安装残留", user_profile / "Downloads", "glob"),
    ]


def normalize_path(raw: str | Path) -> Path:
    return Path(raw).expanduser().resolve()


def is_protected_path(path: Path) -> bool:
    try:
        resolved = path.resolve()
    except OSError:
        return True
    return any(resolved == prefix or prefix in resolved.parents for prefix in PROTECTED_PREFIXES)


def is_excluded(path: Path, excludes: Iterable[Path]) -> bool:
    try:
        resolved = path.resolve()
    except OSError:
        return True
    for exclude in excludes:
        if resolved == exclude or exclude in resolved.parents:
            return True
    return False


def iter_tree_files(target: ScanTarget) -> Iterable[Path]:
    if not target.path.exists():
        logging.info("跳过不存在的目标: %s", target.path)
        return []
    if target.name == "Firefox 缓存":
        files: list[Path] = []
        for profile_dir in target.path.glob("*.default*"):
            cache_dir = profile_dir / "cache2"
            if cache_dir.exists():
                files.extend(p for p in cache_dir.rglob("*") if p.is_file())
        return files
    return (p for p in target.path.rglob("*") if p.is_file())


def iter_glob_files(target: ScanTarget) -> Iterable[Path]:
    if not target.path.exists():
        logging.info("跳过不存在的目标: %s", target.path)
        return []
    files: list[Path] = []
    for pattern in GLOB_PATTERNS:
        files.extend(p for p in target.path.rglob(pattern) if p.is_file())
    return files


def collect_candidates(targets: Iterable[ScanTarget], excludes: list[Path]) -> list[FileCandidate]:
    candidates: list[FileCandidate] = []
    seen: set[Path] = set()

    for target in targets:
        iterator = iter_tree_files(target) if target.mode == "tree" else iter_glob_files(target)
        for path in iterator:
            try:
                resolved = path.resolve()
                if resolved in seen:
                    continue
                if is_excluded(resolved, excludes):
                    continue
                if is_protected_path(resolved) and target.name not in {"Windows 更新缓存", "回收站", "Windows Temp"}:
                    continue
                size = path.stat().st_size
                candidates.append(FileCandidate(path=resolved, size=size, source=target.name))
                seen.add(resolved)
            except (FileNotFoundError, PermissionError, OSError) as exc:
                logging.warning("跳过文件失败: %s (%s)", path, exc)
    return candidates


def remove_readonly(func, path, _exc_info) -> None:
    os.chmod(path, stat.S_IWRITE)
    func(path)


def delete_file(candidate: FileCandidate) -> tuple[bool, str | None]:
    try:
        candidate.path.unlink()
        return True, None
    except IsADirectoryError:
        try:
            shutil.rmtree(candidate.path, onerror=remove_readonly)
            return True, None
        except (PermissionError, OSError) as exc:
            return False, str(exc)
    except (PermissionError, OSError) as exc:
        return False, str(exc)


def format_size(num_bytes: int) -> str:
    size = float(num_bytes)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024 or unit == "TB":
            return f"{size:.2f} {unit}"
        size /= 1024
    return f"{num_bytes} B"


def print_summary(candidates: list[FileCandidate], clean: bool) -> None:
    total_size = sum(item.size for item in candidates)
    action = "将删除" if not clean else "准备删除"
    logging.info("%s %d 个文件，预计释放 %s", action, len(candidates), format_size(total_size))
    for item in candidates:
        logging.info("[%s] %s (%s)", item.source, item.path, format_size(item.size))


def execute_cleanup(candidates: list[FileCandidate]) -> None:
    deleted_count = 0
    deleted_size = 0
    for item in candidates:
        success, error = delete_file(item)
        if success:
            deleted_count += 1
            deleted_size += item.size
            logging.info("已删除: %s", item.path)
        else:
            logging.error("删除失败: %s (%s)", item.path, error)
    logging.info("删除完成: %d 个文件，释放 %s", deleted_count, format_size(deleted_size))


def main() -> int:
    args = parse_args()
    log_path = setup_logging(args.log_file)
    excludes = [normalize_path(item) for item in args.exclude]
    targets = build_targets()

    logging.info("运行模式: %s", "CLEAN" if args.clean else "DRY-RUN")
    logging.info("日志文件: %s", log_path)
    if excludes:
        for item in excludes:
            logging.info("排除路径: %s", item)

    candidates = collect_candidates(targets, excludes)
    print_summary(candidates, args.clean)

    if not args.clean:
        logging.info("当前为 dry-run，未执行删除。使用 --clean 或 --force 才会实际删除。")
        return 0

    execute_cleanup(candidates)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        logging.error("用户中断执行。")
        raise SystemExit(130)
