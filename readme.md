# Segdesign

## 描述
本程序用于对蛋白质单链进行设计，在保留原蛋白质功能的基础上，为其添加新的特性。

## 主要模块
### 1. hmmer模块
基于hmmer工具对输入的蛋白质进行同源搜索，通过同源性分析以及蛋白质二级结构情况最终给出推荐的可设计区域。

### 2. rfdiffusion模块
使用rfdiffusion对蛋白质的设计区域进行定向的三维结构重塑，通过设置的二级结构比例阈值，筛选出符合条件的蛋白质骨架。

### 3. MPNN模块
使用ProteinMPNN（ThermoMPNN或SoluMPNN等）工具对rfdiffusion生成的蛋白质结构的设计区域进行序列预测，根据global_score数值按比例筛选序列，然后使用mmseqs2进行聚类分析，挑选出每个簇的代表序列。

### 4. esmfold模块
使用ESMFold对蛋白质序列进行三维结构预测，可以指定pLDDT筛选阈值，筛选符合条件的序列，最后生成综合性报告。

## 安装说明
### 依赖前提
提前安装好 **anaconda** 和 **dos2unix**。

### 环境构建脚本
各模块的环境构建方法已封装为shell脚本，位于项目根目录，包括：
- hmmer_env.sh
- rfdiffusion_env.sh
- MPNN_env.sh
- esmfold_env.sh
- esmfold_report_env.sh

### 安装步骤
1. 编辑所有shell脚本，将其中的 `CONDA_PATH` 变量改为你的anaconda安装路径。
2. 以hmmer_env.sh为例，运行以下命令构建环境（其他脚本运行方式相同）：
```shell
dos2unix hmmer_env.sh
chmod +x hmmer_env.sh
./hmmer_env.sh
3. 若安装遇到问题，请根据自身电脑配置自行调整。

####补充说明（dssp 4.5.5 问题修复）
dssp 4.5.5 存在兼容问题，涉及环境（hmmer、SE3nv、esmfold_report）需按以下步骤处理：
- 运行 conda env list 查看虚拟环境所在路径。
- 进入 <环境名称>/share/libcifpp 目录。
- 找到 omponents.cif.gz 压缩包并解压。
- 将解压后的文件命名为 omponents.cif，保留在 libcifpp 目录下。
- 删除原压缩包 omponents.cif.gz。