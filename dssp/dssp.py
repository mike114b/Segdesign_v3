import argparse
import os
import subprocess
#import shlex

import sys
from pathlib import Path

# 1. 获取当前脚本(dssp.py)的绝对路径
current_file = Path(__file__).absolute()
# 2. 获取项目根目录（Segdesign_v2）：dssp.py → dssp/ → Segdesign_v2/
project_root = current_file.parent.parent
# 3. 将项目根目录加入Python搜索路径（最关键！）
sys.path.insert(0, str(project_root))  # insert(0)确保优先搜索根目录

# ########## 原导入语句（修改后） ##########
# 此时Python能识别dssp/为包，绝对导入生效
from dssp.dsspcsv import dssp_to_csv

def parse_args():
    parser = argparse.ArgumentParser(description='DSSP analysis',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input_path', '-i', type=str, help='The path of the pdb file or the folder containing the pdb file')
    parser.add_argument('--output_folder', '-o', type=str, help='The folder where the output file is stored')

    return parser.parse_args()

def run_dssp(input_path,output_path):
    output_path_dir = output_path.rsplit('/', 1)[0]
    if not os.path.exists(output_path_dir):
        os.makedirs(output_path_dir,exist_ok=True)

    cmd = ['mkdssp', input_path, output_path]

    # Print command for verification
    #print("=" * 60)
    print("Executing dssp")
    #print(' '.join(shlex.quote(arg) for arg in cmd))
    #print("=" * 60)
    #print()
    # Run command
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=False,  # Show real-time output
            text=True
        )
        print("\n dssp completed successfully!")

    except subprocess.CalledProcessError as e:
        print(f"\n RFdiffusion failed with return code {e.returncode}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print("\n Error: 'python' command or script not found. Check your environment.", file=sys.stderr)
        sys.exit(1)
    return


def main():
    args = parse_args()
    input_path = os.path.expanduser(args.input_path)
    output_path = os.path.expanduser(args.output_folder)
    dssp_folder = os.path.join(output_path,'dssp_files')
    csv_folder = os.path.join(output_path,'csv_files')
    if not os.path.exists(dssp_folder):  ##新建文件夹
        os.makedirs(dssp_folder, exist_ok=True)
    if not os.path.exists(csv_folder):
        os.makedirs(csv_folder, exist_ok=True)
    if os.path.isdir(input_path):
        print('The input path is a folder')
        filenames = os.listdir(input_path)
        for filename in filenames:
            file_name = os.path.splitext(filename)[0]
            file_path = os.path.join(input_path, filename)
            out_dssp_path = os.path.join(dssp_folder, file_name + '.dssp')
            out_csv_path = os.path.join(csv_folder, file_name + '.csv')
            run_dssp(file_path, out_dssp_path)
            dssp_to_csv(out_dssp_path, out_csv_path)

    else:
        print('The input path is a file')
        filename = os.path.basename(input_path)
        file_name = os.path.splitext(filename)[0]
        out_dssp_path = os.path.join(dssp_folder, file_name + '.dssp')
        out_csv_path = os.path.join(csv_folder, file_name + '.csv')
        run_dssp(input_path, out_dssp_path)
        dssp_to_csv(out_dssp_path, out_csv_path)

    print('\nDSSP module has finished running\n')
    return

if __name__ == '__main__':
    main()



