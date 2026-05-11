#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Directory Tree Generator
目录树生成工具 - 递归遍历目录结构并生成JSON格式的树形结构

功能:
- 接受用户输入的绝对路径
- 递归遍历所有文件夹和文件
- 生成JSON格式的树形结构
- 将结果保存到文件并输出到控制台
"""

import json
import sys
from pathlib import Path
from datetime import datetime
from typing import Optional


def build_tree(path: Path) -> dict:
    """
    递归构建目录树

    Args:
        path: 要遍历的路径

    Returns:
        包含目录树结构的字典
    """
    node = {
        "name": path.name,
        "type": "directory" if path.is_dir() else "file",
        "path": str(path.absolute())
    }

    # 如果是目录，递归添加子项
    if path.is_dir():
        try:
            children = []
            for item in sorted(path.iterdir()):
                try:
                    children.append(build_tree(item))
                except (PermissionError, OSError) as e:
                    # 记录无法访问的项目，但继续遍历
                    print(f"Warning: Cannot access {item} - {type(e).__name__}", file=sys.stderr)
            node["children"] = children
        except (PermissionError, OSError) as e:
            # 目录权限拒绝
            print(f"Warning: Cannot read directory {path} - {type(e).__name__}", file=sys.stderr)
            node["children"] = []

    return node


def generate_tree_from_path(input_path: str) -> Optional[dict]:
    """
    从输入路径生成目录树

    Args:
        input_path: 用户输入的路径字符串

    Returns:
        目录树字典，或路径无效时返回 None
    """
    try:
        path = Path(input_path)

        # 验证路径存在
        if not path.exists():
            print(f"Error: Path does not exist: {input_path}", file=sys.stderr)
            return None

        # 构建树形结构
        print(f"Scanning directory: {path.absolute()}")
        tree = build_tree(path)
        return tree

    except (ValueError, TypeError) as e:
        print(f"Error: Invalid path format - {e}", file=sys.stderr)
        return None


def save_to_json(tree: dict, output_dir: Path = None) -> bool:
    """
    将树形结构保存为JSON文件

    Args:
        tree: 目录树字典
        output_dir: 输出目录，默认为当前脚本所在目录

    Returns:
        是否保存成功
    """
    try:
        if output_dir is None:
            output_dir = Path(__file__).parent

        # 确保输出目录存在
        output_dir.mkdir(parents=True, exist_ok=True)

        # 生成带时间戳的文件名
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = output_dir / f"directory_tree_{timestamp}.json"

        # 保存为JSON，使用缩进和中文编码支持
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(tree, f, indent=2, ensure_ascii=False)

        print(f"\n[OK] File saved: {output_file}")
        return True

    except (IOError, OSError) as e:
        print(f"Error: Unable to save file - {e}", file=sys.stderr)
        return False


def print_tree_to_console(tree: dict) -> None:
    """
    将树形结构输出到控制台

    Args:
        tree: 目录树字典
    """
    print("\n" + "="*60)
    print("Generated directory tree (JSON format):")
    print("="*60)
    print(json.dumps(tree, indent=2, ensure_ascii=False))
    print("="*60 + "\n")


def main() -> None:
    """
    主程序 - 交互式界面
    """
    print("="*60)
    print("Directory Tree Generator")
    print("="*60)
    print("Feature: Recursively traverse all files and folders in a path")
    print("Output: Directory tree in JSON format")
    print("-"*60 + "\n")

    while True:
        # 获取用户输入
        user_input = input("Please enter the absolute path (type 'exit' to quit): ").strip()

        # 检查退出命令
        if user_input.lower() in ('exit', 'quit', 'q'):
            print("Program exited")
            break

        # 验证输入是否为空
        if not user_input:
            print("Error: Path cannot be empty, please try again\n")
            continue

        # 生成目录树
        tree = generate_tree_from_path(user_input)

        if tree is None:
            print("Failed to generate directory tree, please check if the path is correct\n")
            continue

        # 输出到控制台
        print_tree_to_console(tree)

        # 保存到文件
        output_dir = Path(__file__).parent
        if save_to_json(tree, output_dir):
            print("[OK] Operation completed")
        else:
            print("[FAIL] Save failed, but the tree structure is displayed above")

        print()


if __name__ == "__main__":
    main()
