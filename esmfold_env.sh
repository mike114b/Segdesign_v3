#!/bin/bash

set -e  # 出错立即退出，避免静默失败
# ESMFold CUDA 12.4 分步安装脚本

echo "=== ESMFold 安装脚本 for CUDA 12.4 ==="

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
echo "创建基础环境..."
echo "create -n esmfold python=3.9 -y"
conda create -n esmfold python=3.9 -y
echo "进入虚拟环境..."
echo "conda activate esmfold"
conda activate esmfold

#2.安装前置
echo "安装前置"
echo "pip install setuptools==59.5.0 --no-cache-dir"
pip install setuptools==59.5.0 --no-cache-dir
echo "conda install conda-forge::openmm=7.5.1 -y"
conda install conda-forge::openmm=7.5.1 -y
echo "conda install -c conda-forge pdbfixer einops fairscale omegaconf=2.1.1 hydra-core=1.1.0 pytest -y "
conda install -c conda-forge pdbfixer einops fairscale omegaconf=2.1.1 hydra-core=1.1.0 pytest -y 

# 3. 安装生物信息学工具
echo "安装生物信息学工具..."
echo "conda install -c conda-forge -c bioconda hhsuite=3.3.0 -y"
conda install -c conda-forge -c bioconda hhsuite=3.3.0 -y
echo "conda install bioconda::hmmer==3.3.2 -y"
conda install bioconda::hmmer==3.3.2 -y
echo "conda install bioconda::kalign2==2.04 -y"
conda install bioconda::kalign2==2.04 -y

# 4. 安装 PyTorch 1.12 (CUDA 11.8版本，但兼容 CUDA 12.4)
echo "安装 PyTorch 1.12..."
echo "pip install torch==2.4.0 torchvision==0.19.0 torchaudio==2.4.0 --index-url https://download.pytorch.org/whl/cu124"
pip install torch==2.4.0 torchvision==0.19.0 torchaudio==2.4.0 --index-url https://download.pytorch.org/whl/cu124

echo "降级pip为24.0版本"
echo "install pip==24.0"
pip install pip==24.0

# 5. 安装 pip 包
echo "安装 pip 包..."
echo "pip install \
biopython==1.79 \
dm-tree==0.1.6 \
ml-collections==0.1.0 \
PyYAML==5.4.1 \
requests==2.26.0 \
scipy  \
tqdm==4.62.2 \
typing-extensions==4.8.0 \
pytorch_lightning==1.6.5 \
wandb==0.12.21 \
panda \
modelcif \
"
pip install \
biopython==1.79 \
dm-tree==0.1.6 \
ml-collections==0.1.0 \
PyYAML==5.4.1 \
requests==2.26.0 \
scipy  \
tqdm==4.62.2 \
typing-extensions==4.8.0 \
pytorch_lightning==1.6.5 \
wandb==0.12.21 \
panda \
modelcif \


#7.安装esmfold
echo "安装esmfold"
echo "pip install 'fair-esm[esmfold]'"
pip install "fair-esm[esmfold]"

# 8. 安装 NVIDIA dllogger
echo "安装 NVIDIA dllogger..."
echo "pip install git+https://github.com/NVIDIA/dllogger.git"
pip install git+https://github.com/NVIDIA/dllogger.git

#9.安装openfold
echo "安装openfold"
echo "cd openfold"
cd openfold
echo "pip install --no-build-isolation ."
pip install . --no-build-isolation
echo "pip uninstall deepspeed -y"
pip uninstall deepspeed -y
