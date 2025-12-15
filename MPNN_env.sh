#!/bin/bash

set -e  # 出错立即退出，避免静默失败


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

TARGET_DIR="./ProteinMPNN"
if [ ! -d "${TARGET_DIR}" ]; then
    # 不存在则创建
   echo "=== 克隆ProteinMPNN库到本地 ==="
   echo "git clone https://github.com/dauparas/ProteinMPNN.git"
   git clone https://github.com/dauparas/ProteinMPNN.git

else
    echo "ProteinMPNN库已存在，跳过"
fi

echo "=== MPNN模块环境的安装 ==="
echo "创建环境..."
echo "conda create -n mpnn python=3.9 -y"
conda create -n mpnn python=3.9 -y
echo "进入虚拟环境..."
echo "conda activate mpnn"
conda activate mpnn

echo "安装pytorch 2.4.0"
echo "pip install torch==2.4.0 torchvision==0.19.0 torchaudio==2.4.0 --index-url https://download.pytorch.org/whl/cu124"
pip install torch==2.4.0 torchvision==0.19.0 torchaudio==2.4.0 --index-url https://download.pytorch.org/whl/cu124

echo "安装biopython库"
echo "conda install biopython -y"
conda install biopython -y
echo "安装mmseqs2"
echo "conda install -c conda-forge -c bioconda mmseqs2 -y"
conda install -c conda-forge -c bioconda mmseqs2 -y