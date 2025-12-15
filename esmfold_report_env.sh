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

echo "=== esmfold_report模块环境的安装 ==="
echo "创建环境..."
echo "conda create -n esmfold_report python=3.9 -y"
conda create -n esmfold_report python=3.9 -y
echo "conda activate esmfold_report"
conda activate esmfold_report
echo "conda install biopython -y"
conda install biopython -y
echo "conda install -c conda-forge dssp=4.5.5 -y"
conda install -c conda-forge dssp=4.5.5 -y
echo "conda install -c conda-forge pandas"
conda install -c conda-forge pandas
echo "conda install conda-forge::biotite"
conda install conda-forge::biotite