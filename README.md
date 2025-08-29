# Fruid / 基于rust+julia的流体模拟引擎

该引擎以安全而具有良好性能的rust作为后端、Julia的[Makie](https://github.com/MakieOrg/Makie.jl)库作为前端，旨在实现一个具有基本功能的、性能高、易于交互的简单流体模拟引擎。  
本模拟引擎的一部分实现思路参考自[该项目](https://github.com/Nicholas-L-Johnson/flip-card)。  
*注：这只是个开发练习项目，并不具备专业水准。*

## 实现原理
项目基于FLIP(Fluid-Implicit-Particle)方法的流体模拟系统其中重要组成部分及原理如下：

# `lib.rs`：

1. 核心数据结构：
- `FlipFluid` 结构体：包含了流体模拟所需的所有数据，如粒子位置、速度、密度等
- `CellType` 枚举：定义了三种单元格类型：空气(AIR_CELL)、流体(FLUID_CELL)和固体(SOLID_CELL)
- `Scene` 结构体：管理整个模拟场景，包括重力、时间步长等参数

2. 主要模拟原理：
- 使用混合欧拉-拉格朗日方法，结合了网格法和粒子法的优点
- 粒子用于追踪流体，网格用于求解压力和非压缩性约束
- 模拟过程包括：
  - 粒子积分（更新位置和速度）
  - 粒子碰撞处理
  - 速度场在网格和粒子之间的转换
  - 求解不可压缩性约束
  - 密度更新

3. 关键算法：
- `integrateParticles`：更新粒子位置和速度
- `pushParticlesApart`：处理粒子间的排斥力，避免过度聚集
- `solveIncompressibility`：求解压力方程，确保流体不可压缩性
- `transferVelocities`：在网格和粒子之间转换速度场
- `handleParticleCollisions`：处理边界碰撞

4. C语言接口：
- 提供了一系列带有`#[unsafe(no_mangle)]`标记的函数
- 允许从C语言调用Rust实现的流体模拟功能
- 包括场景创建、销毁、模拟步进等功能

FLIP方法的优点：
- 比纯粒子法计算效率更高
- 比纯网格法能更好地保持细节
- 数值耗散较小，可以很好地保持流体的动态特性

# `main.jl`

- 自动检查和安装必要的Julia包
- 主要使用Makie进行可视化，TOML进行配置文件解析
- 使用`ccall`函数调用编译好的Rust动态库
- 包装了几个主要的接口：
  - setup_scene：初始化模拟场景
  - destroy_scene：清理场景
  - simulate_step：进行一步模拟
  - set_gravity：设置重力
  - get_output：获取模拟结果
- 使用Makie.jl创建可视化界面
- 使用heatmap热力图方式显示流体分布
- Observable用于实现实时更新
- 使用record函数录制模拟过程
- 每一帧都：
  - 调用Rust进行模拟计算
  - 获取新的模拟数据
  - 更新可视化显示
- 最终生成MP4格式的视频文件

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
- 实现同三维模型的交互
