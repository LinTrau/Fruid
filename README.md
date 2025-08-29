# Fruid / 基于rust+julia的流体模拟引擎

该引擎以安全而具有良好性能的rust作为后端、Julia的[Makie](https://github.com/MakieOrg/Makie.jl)库作为前端，旨在实现一个具有基本功能的、性能高、易于交互的简单流体模拟引擎。  
本模拟器的一部分实现思路参考自[该项目](https://github.com/Nicholas-L-Johnson/flip-card)。  
*注：这只是个开发练习项目，并不具备专业水准。*

## **注意**
您需要安装rust和julia才能使用。

## 使用方法
```bash
#克隆仓库：
git clone https://github.com/LinTrau/Fruid

#进入根目录
cd Fruid

#构建模拟库
cargo build --release

#运行脚本
julia julia/main.jl
```
请在config.toml中输入模拟参数，模拟视频结果将生成于项目根目录的output文件夹。  
注意，如果您是**windows**平台，请将`julia/main.jl`中`const LIB_PATH = abspath(joinpath(SCRIPT_DIR, "../target/release/libfluid_sim.so"))`路径的`.so`改为`.dll`文件。  

## 开发目标
- 实现更多模拟算法
- 实现一个能动态调整参数和设置保存结果的界面
- 加入更多可模拟物理参数
- 提供初步的分析数据

