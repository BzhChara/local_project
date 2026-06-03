import os
import io
import pandas as pd
from pywinauto import Desktop
from PIL import Image, ImageOps, ImageEnhance
import ddddocr
from openpyxl.styles import Font

# 1. 初始化
det = ddddocr.DdddOcr(beta=False, show_ad=False)

# --- 坐标参数 ---
Y_STARTS = [135, 280]
X_START = 98
WIDTH = 52
HEIGHT = 54
X_GAP = 54

def ocr_single_line(img_line, debug_path=None):
    # 1. 左右各切掉 2 像素，避开深色竖线干扰
    w, h = img_line.size
    img_line = img_line.crop((2, 2, w - 2, h-2)) 
    
    # 2. 放大 5 倍并转灰度
    img_line = img_line.convert('L').resize((img_line.width * 5, img_line.height * 5), Image.LANCZOS)
    
    # 3. 锐化
    img_line = ImageEnhance.Sharpness(img_line).enhance(2.0)
    
    # 4. 二值化
    img_line = img_line.point(lambda x: 0 if x < 155 else 255)
    
    # 5. 保存 OCR 实际看到的处理后的图
    if debug_path:
        img_line.save(debug_path)
    
    img_byte_arr = io.BytesIO()
    img_line.save(img_byte_arr, format='PNG')
    return det.classification(img_byte_arr.getvalue())

def run_task_with_processed_images():
    debug_folder = "processed_verify_crops"
    if not os.path.exists(debug_folder):
        os.makedirs(debug_folder)
        
    try:
        app = Desktop(backend="win32")
        win = app.window(title_re=".*softwinner.*")
        win.set_focus()
        full_img = win.capture_as_image()
        print(f"正在识别并将处理后的【黑白稿】保存至 {debug_folder}...")

        all_results = []
        line_h = HEIGHT // 3  
        labels = ["原始值", "校准值", "基准值"]
        suffixes = ["_row1_raw", "_row2_cal", "_row3_base"]

        for row_idx, row_y in enumerate(Y_STARTS):
            for i in range(12):
                ch_num = row_idx * 12 + i + 1
                x = X_START + i * X_GAP
                roi = full_img.crop((x, row_y, x + WIDTH, row_y + HEIGHT))
                row_data = {"通道号": f"CH{ch_num:02d}"}
                
                for line_idx in range(3):
                    top = line_idx * line_h
                    bottom = (line_idx + 1) * line_h
                    line_img = roi.crop((0, top, WIDTH, bottom))
                    
                    # 生成带路径的文件名
                    img_filename = f"CH{ch_num:02d}{suffixes[line_idx]}.png"
                    img_debug_path = os.path.join(debug_folder, img_filename)
                    
                    # 传入路径，让识别函数保存处理后的图
                    text = ocr_single_line(line_img, debug_path=img_debug_path)
                    
                    clean_text = "".join([c for c in text if c.isdigit() or c == '.'])
                    row_data[labels[line_idx]] = clean_text

                print(f"通道 {ch_num:02d} -> 1:{row_data['原始值']:6} | 2:{row_data['校准值']:6} | 3:{row_data['基准值']:6}")
                all_results.append(row_data)

        # 导出 Excel
        df = pd.DataFrame(all_results)
        file_name = "Sigma_黑白稿对比结果.xlsx"
        with pd.ExcelWriter(file_name, engine='openpyxl') as writer:
            df.to_excel(writer, index=False)
            ws = writer.sheets['Sheet1']
            for row in ws.iter_rows():
                for cell in row:
                    cell.font = Font(name='微软雅黑', size=11)

        print(f"\n完成！请对照 '{debug_folder}' 文件夹里的黑白图片检查。")

    except Exception as e:
        print(f"报错: {e}")

if __name__ == "__main__":
    run_task_with_processed_images()