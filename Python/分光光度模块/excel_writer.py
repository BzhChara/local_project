import csv
import os
import time

from openpyxl import Workbook
from openpyxl.styles import Font

from app_paths import get_app_dir


EXCEL_FONT_NAME = "微软雅黑"
TEMP_CSV_DIR_NAME = "_temp_csv"


def check_output_available(path: str):
    if not os.path.exists(path):
        return True

    temp_output = path + ".lockcheck"
    try:
        os.rename(path, temp_output)
        os.rename(temp_output, path)
    except OSError:
        print(f"输出文件被占用，请先关闭: {os.path.abspath(path)}")
        return False

    return True


def create_run_output_dir():
    run_time = time.strftime("%H%M%S")
    output_dir = os.path.join(get_app_dir(), "data", f"record_{run_time}")
    os.makedirs(output_dir, exist_ok=True)
    return output_dir


def resolve_output_path(output_dir: str, file_name: str):
    return os.path.join(output_dir, os.path.basename(file_name))


def raw_csv_path(raw_excel_path: str):
    temp_dir = os.path.join(os.path.dirname(raw_excel_path), TEMP_CSV_DIR_NAME)
    csv_name = f"{os.path.splitext(os.path.basename(raw_excel_path))[0]}.csv"
    return os.path.join(temp_dir, csv_name)


def _selected_indices(target_indices=None):
    return list(range(24)) if target_indices is None else list(target_indices)


def create_raw_csv(raw_excel_path: str, target_indices=None):
    selected_indices = _selected_indices(target_indices)
    csv_path = raw_csv_path(raw_excel_path)
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    file = open(csv_path, "w", newline="", encoding="utf-8-sig")
    writer = csv.writer(file)
    writer.writerow(["时间"] + [f"CH{idx + 1}" for idx in selected_indices])
    file.flush()
    return file, writer, csv_path


def append_raw_csv_row(file, writer, values: list[float], target_indices=None):
    selected_indices = _selected_indices(target_indices)
    row = [time.strftime("%H:%M:%S")] + [f"{values[idx]:.6f}" for idx in selected_indices]
    writer.writerow(row)
    file.flush()


def create_raw_excel(path: str):
    wb = Workbook()
    ws = wb.active
    ws.title = "Sheet1"

    headers = ["时间"] + [f"CH{i}" for i in range(1, 25)]
    ws.append(headers)
    ws.freeze_panes = "B2"

    for cell in ws[1]:
        cell.font = Font(name=EXCEL_FONT_NAME, size=12, bold=True)

    ws.column_dimensions["A"].width = 12
    for col_idx in range(2, 26):
        col_letter = ws.cell(row=1, column=col_idx).column_letter
        ws.column_dimensions[col_letter].width = 12

    wb.save(path)
    return wb, ws


def raw_csv_to_excel(csv_path: str, excel_path: str):
    wb = Workbook()
    ws = wb.active
    ws.title = "Sheet1"
    ws.freeze_panes = "B2"

    with open(csv_path, "r", newline="", encoding="utf-8-sig") as file:
        reader = csv.reader(file)
        for row_index, row in enumerate(reader, start=1):
            if row_index == 1:
                ws.append(row)
            else:
                values = [row[0]] + [float(value) if value else None for value in row[1:]]
                ws.append(values)

            for cell in ws[row_index]:
                cell.font = Font(name=EXCEL_FONT_NAME, size=12, bold=row_index == 1)
            for cell in ws[row_index][1:]:
                cell.number_format = "0.000000"

    ws.column_dimensions["A"].width = 12
    for col_idx in range(2, ws.max_column + 1):
        col_letter = ws.cell(row=1, column=col_idx).column_letter
        ws.column_dimensions[col_letter].width = 12

    wb.save(excel_path)
    return wb, ws


def create_stable_excel(path: str, target_indices=None):
    wb = Workbook()
    ws = wb.active
    ws.title = "Sheet1"
    refresh_stable_excel(wb, ws, path, [[] for _ in range(24)], target_indices)
    return wb, ws


def _stable_history(stable_values: list, idx: int):
    value = stable_values[idx]
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def _stats(values: list[float]):
    if not values:
        return [None, None, None, None]
    avg_value = sum(values) / len(values)
    std_value = (sum((value - avg_value) ** 2 for value in values) / len(values)) ** 0.5
    cv_value = std_value / avg_value if avg_value else 0
    return [min(values), max(values), avg_value, cv_value]


def _write_cell(ws, row: int, column: int, value, bold: bool = False, cv: bool = False):
    cell = ws.cell(row=row, column=column, value=value)
    cell.font = Font(name=EXCEL_FONT_NAME, size=12, bold=bold)
    if isinstance(value, (int, float)):
        cell.number_format = "0.00%" if cv else "0.000000"
    return cell


def _reset_sheet(ws):
    if ws.max_row:
        ws.delete_rows(1, ws.max_row)
    if ws.max_column:
        ws.delete_cols(1, ws.max_column)


def refresh_stable_excel(wb: Workbook, ws, path: str, stable_values: list, target_indices=None):
    selected_indices = _selected_indices(target_indices)
    histories = {idx: _stable_history(stable_values, idx) for idx in selected_indices}
    max_count = max((len(values) for values in histories.values()), default=0)
    measurement_count = max(max_count, 1)
    stats_col = measurement_count + 3
    bottom_stats_start = len(selected_indices) + 3

    _reset_sheet(ws)

    _write_cell(ws, 1, 1, "通道", bold=True)
    for attempt in range(measurement_count):
        _write_cell(ws, 1, 2 + attempt, f"第{attempt + 1}次", bold=True)

    for offset, label in enumerate(["最小值", "最大值", "平均数", "CV"]):
        _write_cell(ws, 1, stats_col + offset, label, bold=True)

    for row_idx, idx in enumerate(selected_indices, start=2):
        history = histories[idx]
        _write_cell(ws, row_idx, 1, f"CH{idx + 1}", bold=True)
        for attempt, value in enumerate(history, start=1):
            _write_cell(ws, row_idx, 1 + attempt, value)

        for offset, value in enumerate(_stats(history)):
            _write_cell(ws, row_idx, stats_col + offset, value, cv=offset == 3)

    for offset, label in enumerate(["最小值", "最大值", "平均数", "CV"]):
        _write_cell(ws, bottom_stats_start + offset, 1, label, bold=True)
        for attempt in range(measurement_count):
            values = [history[attempt] for history in histories.values() if len(history) > attempt]
            _write_cell(ws, bottom_stats_start + offset, 2 + attempt, _stats(values)[offset], cv=offset == 3)

    ws.column_dimensions["A"].width = 12
    for col_idx in range(2, stats_col + 4):
        col_letter = ws.cell(row=1, column=col_idx).column_letter
        ws.column_dimensions[col_letter].width = 14

    ws.column_dimensions["A"].width = 12

    wb.save(path)
