import subprocess
import os
import logging
from typing import Dict, Optional, List
import shlex
import argparse
from pathlib import Path
import yaml
import sys
import threading

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(stream=sys.stdout),
        logging.FileHandler('module_runner.log', encoding='utf-8')
    ]
)
logger = logging.getLogger(__name__)

# 配置项（可根据实际情况修改）
CONFIG = {
    "MODULES":{
        'hmmer': {"DIR":'./hmmer'},
        'rf_diffusion': {"DIR":'./rfdiffusion'},
        'MPNN': {"DIR":'./MPNN'},
        'esmfold': {"DIR":'./esmfold'},
        'esmfold_report': {"DIR":'./esmfold'},
        'dssp': {"DIR":'./dssp'},
        'cluster_analysis':{"DIR":'./MPNN'},
    },

}



class ModuleRunnerError(Exception):
    """模块运行器自定义异常"""
    pass


def validate_environment(env_name: str) -> bool:
    """验证Conda环境是否存在"""
    conda_info_cmd = [
        f"{CONFIG['MINICONDA_PATH']}/bin/conda",
        "info",
        "--envs"
    ]

    try:
        result = subprocess.run(
            conda_info_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
            timeout=30
        )
        # 检查环境是否在输出中（支持完整名称匹配）
        return any(f"*{env_name}" in line or f"  {env_name} " in line for line in result.stdout.splitlines())
    except subprocess.TimeoutExpired:
        logger.warning(f"验证环境 {env_name} 超时")
        return False
    except subprocess.CalledProcessError as e:
        logger.error(f"验证环境失败: {e.stderr}")
        return False


def validate_module(module_name: str) -> str:
    """验证模块是否存在并返回完整路径"""
    if module_name not in CONFIG['MODULES']:
        raise ModuleRunnerError(f"模块 {module_name} 未在配置中定义，可用模块: {list(CONFIG['MODULES'].keys())}")

    module_path = os.path.abspath(os.path.join(CONFIG['MODULES'][module_name]['DIR'], f"{module_name}.py"))
    if not os.path.exists(module_path):
        raise ModuleRunnerError(f"模块文件不存在: {module_path}")

    if not os.access(module_path, os.R_OK):
        raise ModuleRunnerError(f"模块文件无读取权限: {module_path}")

    return module_path


def build_command(module_name: str, module_path: str, anaconda_path: str, env_name: str, custom_args: List[str]) -> str:
    """构建安全的执行命令"""


    # 合并默认参数和自定义参数（自定义参数优先级更高）
    #default_args = MODULE_CONFIG[module_name]["default_args"]
    #final_args = default_args + custom_args

    # 安全转义所有参数，防止命令注入
    escaped_args = [shlex.quote(arg) for arg in custom_args]
    args_str = " ".join(escaped_args)

    # 构建命令（使用set -e确保任一命令失败即退出）
    command = f"""
    #!/bin/bash
    set -euo pipefail
    PS1="${{PS1:-}}"
    # 加载conda环境
    if [ -f "{shlex.quote(anaconda_path)}/etc/profile.d/conda.sh" ]; then
        source "{shlex.quote(anaconda_path)}/etc/profile.d/conda.sh"
    elif [ -f "{shlex.quote(anaconda_path)}/bin/activate" ]; then
        source "{shlex.quote(anaconda_path)}/bin/activate"
    else
        echo "找不到conda激活脚本" >&2
        exit 1
    fi

    # 激活环境并运行模块
    conda activate {shlex.quote(env_name)}
    python {shlex.quote(module_path)} {args_str}
    """

    return command
def run_command(command):
    # 创建子进程，捕获标准输出和错误
    print('*'*10)
    print(f"Now starting to execute the command:\n{command}")
    print('*'*10)
    process = subprocess.Popen(
            command,
            shell=True,
            executable="/bin/bash",
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
    # 实时打印输出的函数
    def print_output():
        for line in iter(process.stdout.readline, ''):
            # 移除行尾换行符后打印
            print(line, end='')
            sys.stdout.flush()  # 确保立即显示
        process.stdout.close()
    # 启动输出打印线程
    output_thread = threading.Thread(target=print_output)
    output_thread.daemon = True  # 主程序退出时自动结束线程
    output_thread.start()
    # 等待进程结束
    process.wait()
    # 检查退出状态
    if process.returncode != 0:
        raise RuntimeError(f"Command execution failed，exit code: {process.returncode}")
    return


def run_module(
        module_name: str,
        anaconda_path,
        params,
        retry_count: int = 0
) :
    """
    在指定Conda环境中运行模块（支持重试）

    Args:
        module_name: 模块名称
        args: 模块的命令行参数
        retry_count: 当前重试次数

    Returns:
        退出代码（0表示成功）

    Raises:
        ModuleRunnerError: 模块验证或运行失败时抛出
    """
    # 验证模块
    try:
        module_path = validate_module(module_name)
    except ModuleRunnerError as e:
        logger.error(f"模块验证失败: {e}")
        raise

    # 获取环境名称
    env_name = params['env_name']
    logger.info(f"🚀 启动模块: {module_name} (环境: {env_name}, 路径: {module_path})")

    args = [elem for k, v in params['args'].items() for elem in (f'--{k}', str(v))]
    # 构建命令
    command = build_command(
        module_name=module_name,
        module_path=module_path,
        anaconda_path=os.path.expanduser(anaconda_path),
        env_name=env_name,
        custom_args=list(args)
    )

    run_command(command)
    return



def run_module_old(
        module_name: str,
        anaconda_path,
        params,
        retry_count: int = 0
) -> int:
    """
    在指定Conda环境中运行模块（支持重试）

    Args:
        module_name: 模块名称
        args: 模块的命令行参数
        retry_count: 当前重试次数

    Returns:
        退出代码（0表示成功）

    Raises:
        ModuleRunnerError: 模块验证或运行失败时抛出
    """
    # 验证模块
    try:
        module_path = validate_module(module_name)
    except ModuleRunnerError as e:
        logger.error(f"模块验证失败: {e}")
        raise

    # 获取环境名称
    env_name = params['env_name']
    logger.info(f"🚀 启动模块: {module_name} (环境: {env_name}, 路径: {module_path})")

    args = [elem for k, v in params['args'].items() for elem in (f'--{k}', str(v))]
    # 构建命令
    command = build_command(
        module_name=module_name,
        module_path=module_path,
        anaconda_path=os.path.expanduser(anaconda_path),
        env_name=env_name,
        custom_args=list(args)
    )

    try:
        # 执行命令
        result = subprocess.run(
            command,
            shell=True,
            executable="/bin/bash",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=CONFIG["COMMAND_TIMEOUT"]
        )

        # 记录输出
        logger.info(f"=== 模块 {module_name} 输出 ===")
        if result.stdout:
            logger.info(result.stdout)
        if result.stderr:
            logger.error(f"模块 {module_name} 错误输出: {result.stderr}")

        logger.info(f"模块 {module_name} 退出代码: {result.returncode}")

        # 重试逻辑
        #if result.returncode != 0 and retry_count < CONFIG["MAX_RETRIES"]:
            #retry_count += 1
            #logger.warning(f"模块 {module_name} 运行失败，将进行第 {retry_count}/{CONFIG['MAX_RETRIES']} 次重试...")
            #return run_module(module_name, *args, retry_count=retry_count)

        return result.returncode

    except subprocess.TimeoutExpired:
        error_msg = f"模块 {module_name} 运行超时（{CONFIG['COMMAND_TIMEOUT']}秒）"
        logger.error(error_msg)
        raise ModuleRunnerError(error_msg) from None
    except subprocess.CalledProcessError as e:
        error_msg = f"模块 {module_name} 运行失败: {e.stderr}"
        logger.error(error_msg)
        raise ModuleRunnerError(error_msg) from e
    except Exception as e:
        error_msg = f"模块 {module_name} 运行异常: {str(e)}"
        logger.error(error_msg, exc_info=True)
        raise ModuleRunnerError(error_msg) from e


def read_yaml_file(yaml_path: str) -> dict:
    """
    读取YAML文件并返回字典格式数据

    Args:
        yaml_path: YAML文件的路径（相对路径或绝对路径）

    Returns:
        解析后的字典数据

    Raises:
        FileNotFoundError: 文件不存在
        yaml.YAMLError: YAML格式错误
        PermissionError: 无文件读取权限
    """
    # 转换为Path对象，方便路径处理
    file_path = Path(yaml_path)

    # 检查文件是否存在
    if not file_path.exists():
        raise FileNotFoundError(f"错误：文件不存在 → {yaml_path}")

    # 检查是否是文件（不是目录）
    if not file_path.is_file():
        raise IsADirectoryError(f"错误：{yaml_path} 是目录，不是文件")

    # 读取并解析YAML文件
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            # yaml.safe_load() 避免执行恶意代码，更安全
            data = yaml.safe_load(f)
        return data
    except PermissionError:
        raise PermissionError(f"错误：无权限读取文件 → {yaml_path}")
    except yaml.YAMLError as e:
        raise yaml.YAMLError(f"错误：YAML格式无效 → {e}")
    except Exception as e:
        raise Exception(f"未知错误：{e}")

def global_work_dir_handling(yaml_data):
    work_dir = os.path.expanduser(yaml_data['global parameters']["work_dir"])
    if not os.path.exists(work_dir):
        os.makedirs(work_dir, exist_ok=True)
    if 'hmmer' in yaml_data:
        yaml_data['hmmer']['args']['output_folder'] = os.path.join(work_dir, "hmmer_out")
    if 'rf_diffusion' in yaml_data:
        yaml_data['rf_diffusion']['args']['inference.output_prefix'] = os.path.join(work_dir, "rfdiffusion_out/sample")
    if 'MPNN' in yaml_data:
        yaml_data['MPNN']['args']['output_folder'] = os.path.join(work_dir, "mpnn_out")
    if 'esmfold' in yaml_data:
        yaml_data['esmfold']['args']['output_folder'] = os.path.join(work_dir, "esmfold_out")
    if 'esmfold_report' in yaml_data:
        yaml_data['esmfold_report']['args']['esmfold_folder'] = os.path.join(work_dir, "esmfold_out")
    if 'dssp' in yaml_data:
        yaml_data['dssp']['args']['output_folder'] = os.path.join(work_dir, "dssp_out")
    if 'cluster_analysis' in yaml_data:
        yaml_data['cluster_analysis']['args']['output_folder'] = os.path.join(work_dir, "cluster_analysis_out")
    return yaml_data







if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="通过命令行传入YAML文件路径并读取其内容",
        epilog="示例：python read_yaml_from_cli.py config.yaml"
    )

    # 2. 添加必选参数：yaml文件路径
    parser.add_argument(
        "yaml_file",  # 参数名（位置参数，无需--前缀）
        type=str,
        default="parameter.yaml",
        help="YAML文件的路径（相对路径或绝对路径）"
    )
    args = parser.parse_args()
    try:
        yaml_data = read_yaml_file(args.yaml_file)
        print("✅ YAML文件读取成功！")
        print("📊 解析后的数据：")
        # 格式化输出（可选，更易读）
        print(yaml.dump(yaml_data, allow_unicode=True, sort_keys=False))
    except Exception as e:
        print(f"❌ 读取失败：{e}")
        exit(1)  # 非0退出码表示程序异常

    print(yaml_data)

    if 'work_dir' in yaml_data['global parameters'] and yaml_data['global parameters'] is not None:
        yaml_data = global_work_dir_handling(yaml_data)
    yaml_module = list(yaml_data)
    for module_name in yaml_module:
        if module_name in CONFIG['MODULES']:
            try:
                # 方式1: 运行所有模块
                # results = run_all_modules()

                # 方式2: 单独运行指定模块
                run_module(
                    module_name=module_name,
                    anaconda_path=yaml_data['global parameters']['anaconda_path'],
                    params=yaml_data[module_name]
                )
                # run_module("module2", "--mode", "fast")

                # 方式3: 使用默认参数运行
                # run_module("module2")

            except ModuleRunnerError as e:
                logger.critical(f"程序执行失败: {e}")
                exit(1)
            except KeyboardInterrupt:
                logger.info("程序被用户中断")
                exit(0)
            except Exception as e:
                logger.critical(f"未预期的错误: {str(e)}", exc_info=True)
                exit(1)
        #else:
            #print(f'module{module_name} not enabled')


