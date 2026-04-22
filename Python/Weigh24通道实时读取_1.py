import time
import pandas as pd
from datetime import datetime
from pywinauto import Application

# --- 配置区 ---
EXE_NAME = "Weigh.exe"
OUTPUT_FILE = "24通道原始值.xlsx"
INTERVAL = 3

def is_valid_float(s):
    # 判断是否为我们要的原始值：包含小数点且能转为浮点数
    try:
        s = str(s).strip()
        if "." in s:
            float(s)
            return True
        return False
    except:
        return False

def main():
    try:
        print("正在连接 Weigh.exe...")
        app = Application(backend="uia").connect(path=EXE_NAME)
        dlg = app.top_window()

        all_records = []
        while True:
            current_time = datetime.now().strftime("%H:%M:%S")
            
            # 1. 重新抓取所有控件并排序
            all_controls = dlg.descendants(control_type="Edit") + dlg.descendants(control_type="ComboBox")
            
            # 排除左侧无关区域
            min_left = min([b.rectangle().left for b in all_controls])
            data_controls = [b for b in all_controls if b.rectangle().left > min_left + 150]
            
            # 严格按视觉排序：先从上到下，再从左到右
            data_controls.sort(key=lambda b: (round(b.rectangle().top / 5), b.rectangle().left))

            # 2. 遍历所有控件，只挑出小数
            valid_floats = []
            for ctrl in data_controls:
                try:
                    # 尝试多种方式读取
                    val = ctrl.get_value() or ctrl.window_text() or ""
                    if is_valid_float(val):
                        valid_floats.append(val.strip())
                except:
                    continue

            # 3. 结果对齐
            # 理论上 valid_floats 现在应该正好有 24 个（前12个是第一排，后12个是第二排）
            if len(valid_floats) >= 24:
                final_24 = valid_floats[:24]
                status = "完整采集"
            else:
                # 如果没抓满，用0补齐，防止程序崩溃
                final_24 = valid_floats + ["0"] * (24 - len(valid_floats))
                status = f"数据不足(仅抓到{len(valid_floats)}个)"

            # 4. 保存数据
            row = [current_time] + final_24
            all_records.append(row)
            
            try:
                cols = ['时间'] + [f'CH{i+1}' for i in range(24)]
                df = pd.DataFrame(all_records, columns=cols)
                df.to_excel(OUTPUT_FILE, index=False)
                
                # 打印重点通道核对
                print(f"[{current_time}] {status} | CH1:{final_24[0]} | CH12:{final_24[11]} | CH13:{final_24[12]} | CH24:{final_24[23]}")
            except PermissionError:
                print(f"[{current_time}] Excel 被占用，请关闭文件")

            time.sleep(INTERVAL)

    except Exception as e:
        print(f"发生错误: {e}")

if __name__ == "__main__":
    main()