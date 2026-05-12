from led_debugger.dependency_check import _version_tuple, decode_process_output, parse_requirements


def test_parse_requirements_reads_package_name_and_version_spec(tmp_path) -> None:
    requirements = tmp_path / "requirements.txt"
    requirements.write_text(
        "\n".join(
            [
                "PySide6>=6.7",
                "pyserial>=3.5",
                "",
                "# comment",
                "pytest==8.0",
            ]
        ),
        encoding="utf-8",
    )

    items = parse_requirements(requirements)

    assert [item.name for item in items] == ["PySide6", "pyserial", "pytest"]
    assert items[0].operator == ">="
    assert items[0].required_version == "6.7"
    assert items[2].operator == "=="
    assert items[2].required_version == "8.0"


def test_version_tuple_compares_numeric_parts() -> None:
    assert _version_tuple("6.7.1") > _version_tuple("6.7")
    assert _version_tuple("0.14.0") > _version_tuple("0.13")


def test_decode_process_output_handles_chinese_path_bytes() -> None:
    text = "D:\\Microsoft VS Code Projects\\Python\\上位机程序\\requirements.txt"

    assert decode_process_output(text.encode("gbk")) == text
