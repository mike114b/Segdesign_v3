#!/bin/bash

set -e  # 出错立即退出，避免静默失败
# Hmmer分步安装脚本

echo "=== Hmmer模块环境的安装 ==="

#写入你anaconda的安装路径
CONDA_PATH="/opt/software/anaconda3"
# 加载conda环境
if [ -f "$CONDA_PATH/etc/profile.d/conda.sh" ]; then
    source "$CONDA_PATH/etc/profile.d/conda.sh"
elif [ -f "$CONDA_PATH/bin/activate" ]; then
    source "$CONDA_PATH/bin/activate"
else
    echo "找不到conda激活脚本" >&2
    exit 1
fi
# 1. 创建基础环境
echo "创建环境..."
echo "conda create -n hmmer python=3.9 -y"
conda create -n hmmer python=3.9 -y
echo "进入虚拟环境..."
echo "conda activate hmmer"
conda activate hmmer


# 添加 conda-forge 频道（若未添加）
echo "添加 conda-forge 频道"
conda config --add channels conda-forge
conda config --set channel_priority strict
# 安装 evcouplings（自动安装所有依赖，如 numpy、scipy、pandas 等）

echo "安装hmmer"
echo "conda install bioconda::hmmer==3.3.2 -y"
conda install bioconda::hmmer==3.3.2 -y


#安装pandas库
echo "安装pandas库"
echo "conda install pandas -y"
conda install pandas -y
#安装biopython库
echo "安装biopython库"
echo "conda install biopython -y"
conda install biopython -y
# 安装 hhsuite
echo "安装hhsuite库"
echo "conda install -c conda-forge -c bioconda hhsuite -y"
conda install -c conda-forge -c bioconda hhsuite -y

echo "安装evcouplings库"
echo "pip install evcouplings"
pip install evcouplings

# 安装dssp
echo "安装dssp"
echo "conda install -c conda-forge dssp=4.5.5 -y"
conda install -c conda-forge dssp=4.5.5 -y
