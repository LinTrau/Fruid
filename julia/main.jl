using Pkg

function install_dependencies()
    packages = ["TOML", "Makie", "GLMakie", "Libdl", "FileIO", "Statistics"]
    for pkg in packages
        if !haskey(Pkg.installed(), pkg)
            println("正在安装包: $pkg...")
            Pkg.add(pkg)
            println("包 $pkg 安装完成。")
        end
    end
end

install_dependencies()

using TOML
using Makie
using GLMakie
using Libdl
using FileIO
using Statistics

# 获取 Julia 脚本的目录
const SCRIPT_DIR = @__DIR__

# 构建一个绝对路径来定位 Rust 库文件和配置文件
const LIB_PATH = abspath(joinpath(SCRIPT_DIR,
    "../target/release/libFruid.so"))
const CONFIG_PATH = abspath(joinpath(SCRIPT_DIR, "../config.toml"))

# 检查库文件是否存在
if !isfile(LIB_PATH)
    error("找不到编译后的 Rust 库文件: $(LIB_PATH)。请确保你已成功运行 `cargo build --release`。")
end

const LIB = LIB_PATH

function read_config(config_file::String)
    return TOML.parsefile(config_file)
end

function setup_scene(num_particles::Int32, num_rows::Int32, num_cols::Int32)::Ptr{Cvoid}
    ccall((:setup_scene_wrapper, LIB), Ptr{Cvoid}, (Cint, Cint, Cint), num_particles, num_rows, num_cols)
end

function destroy_scene(scene_ptr::Ptr{Cvoid})
    ccall((:destroy_scene_wrapper, LIB), Cvoid, (Ptr{Cvoid},), scene_ptr)
end

function simulate_step(scene_ptr::Ptr{Cvoid})
    ccall((:simulate_wrapper, LIB), Cvoid, (Ptr{Cvoid},), scene_ptr)
end

function set_gravity(scene_ptr::Ptr{Cvoid}, x::Float32, y::Float32)
    ccall((:set_gravity_wrapper, LIB), Cvoid, (Ptr{Cvoid}, Cfloat, Cfloat), scene_ptr, x, y)
end

function get_output_dims(scene_ptr::Ptr{Cvoid})::Tuple{Int,Int}
    rows = ccall((:get_output_rows, LIB), Cint, (Ptr{Cvoid},), scene_ptr)
    cols = ccall((:get_output_cols, LIB), Cint, (Ptr{Cvoid},), scene_ptr)
    return Int(rows), Int(cols)
end

function get_output(scene_ptr::Ptr{Cvoid})
    rows, cols = get_output_dims(scene_ptr)
    ptr = ccall((:get_output_wrapper, LIB), Ptr{Bool}, (Ptr{Cvoid},), scene_ptr)
    output_array = unsafe_wrap(Array, ptr, (cols, rows), own=true)
    return output_array'
end

function get_density(scene_ptr::Ptr{Cvoid})
    rows, cols = get_output_dims(scene_ptr)
    ptr = ccall((:get_density_wrapper, LIB), Ptr{Cfloat}, (Ptr{Cvoid},), scene_ptr)
    density_array = unsafe_wrap(Array, ptr, (cols, rows), own=false)
    return density_array'
end

function get_u_velocity(scene_ptr::Ptr{Cvoid})
    rows, cols = get_output_dims(scene_ptr)
    ptr = ccall((:get_u_velocity_wrapper, LIB), Ptr{Cfloat}, (Ptr{Cvoid},), scene_ptr)
    u_array = unsafe_wrap(Array, ptr, (cols, rows), own=false)
    return u_array'
end

function get_v_velocity(scene_ptr::Ptr{Cvoid})
    rows, cols = get_output_dims(scene_ptr)
    ptr = ccall((:get_v_velocity_wrapper, LIB), Ptr{Cfloat}, (Ptr{Cvoid},), scene_ptr)
    v_array = unsafe_wrap(Array, ptr, (cols, rows), own=false)
    return v_array'
end

function get_particles(scene_ptr::Ptr{Cvoid})
    num_particles = ccall((:get_num_particles_wrapper, LIB), Cint, (Ptr{Cvoid},), scene_ptr)
    ptr = ccall((:get_particle_positions_wrapper, LIB), Ptr{Cfloat}, (Ptr{Cvoid},), scene_ptr)
    particle_array = unsafe_wrap(Array, ptr, (2, Int(num_particles)), own=false)
    return particle_array
end

function simulation()
    config = read_config(CONFIG_PATH)

    num_particles = config["simulation"]["num_particles"]
    sim_steps = config["simulation"]["sim_steps"]
    gravity = config["simulation"]["gravity"]
    window_width = config["window"]["width"]
    window_height = config["window"]["height"]

    println("从配置文件读取参数：")
    println("粒子数量: ", num_particles)
    println("重力: ", gravity)
    println("模拟步数: ", sim_steps)
    println("窗口尺寸: (", window_width, ", ", window_height, ")")

    output_dir = abspath(joinpath(SCRIPT_DIR, "../output"))
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    gravity_x_str = replace(string(gravity[1]), "." => "p", "-" => "m")
    gravity_y_str = replace(string(gravity[2]), "." => "p", "-" => "m")
    filename = "particles_$(num_particles)_steps_$(sim_steps)_gravity_$(gravity_x_str)_$(gravity_y_str).mp4"
    output_path = joinpath(output_dir, filename)

    scene_ptr = setup_scene(Int32(num_particles), Int32(window_height), Int32(window_width))
    println("场景已初始化。")

    set_gravity(scene_ptr, Float32(gravity[1]), Float32(gravity[2]))

    # 先运行一步模拟以获取初始数据
    simulate_step(scene_ptr)

    # 初始化数据
    density_data = Observable(get_density(scene_ptr))
    particles = Observable(get_particles(scene_ptr))

    # 从粒子数据中提取x和y坐标
    particle_x = @lift($particles[1, :])
    particle_y = @lift($particles[2, :])

    # 设置主题和创建图形
    set_theme!(backgroundcolor=:black)
    fig = Figure(resolution=(800, 600))
    ax = Axis(fig[1, 1], aspect=DataAspect(), title="Fluid Simulation - FLIP Method", xlabel="X", ylabel="Y", xlabelpadding=10, ylabelpadding=10)

    # 设置坐标轴范围
    xlims!(ax, 0, window_width)
    ylims!(ax, 0, window_height)

    # 1. 密度场背景 (使用热图)
    println("密度场数据范围: ", extrema(density_data[]))
    hm = heatmap!(ax, density_data, colormap=:plasma, alpha=0.7, colorrange=(0, maximum(density_data[])))
    Colorbar(fig[1, 2], hm, label="Density")

    # 2. 粒子散点图 (主要的可视化元素)
    println("粒子位置范围: X=", extrema(particle_x[]), " Y=", extrema(particle_y[]))
    scatter!(ax, particle_x, particle_y, color=:cyan, markersize=8, strokewidth=1, strokecolor=:white, alpha=0.8)

    # 添加一些统计信息的文本
    particle_count_text = @lift("Particle numbers: $(length($particle_x))")
    text!(ax, 0.02, 0.98, text=particle_count_text, space=:data, align=(:left, :top), color=:white, fontsize=12)

    display(fig)
    println("开始录制动画...")

    # 录制动画
    record(fig, output_path, 1:sim_steps; framerate=24) do frame
        # 执行模拟步骤
        simulate_step(scene_ptr)

        # 获取新数据
        new_density = get_density(scene_ptr)
        new_particles = get_particles(scene_ptr)

        # 更新Observable数据
        density_data[] = new_density
        particles[] = new_particles

        # 每10帧打印一次进度
        if frame % 10 == 0
            println("帧: $frame / $sim_steps ($(round(frame/sim_steps*100, digits=1))%)")
            println("  最大密度位置: ", argmax(new_density))
            println("  密度范围: $(extrema(new_density))")

            # 检查粒子位置是否有变化
            if frame > 10
                pos_x = new_particles[1, :]
                pos_y = new_particles[2, :]
                println("  X位置范围: $(extrema(pos_x))")
                println("  Y位置范围: $(extrema(pos_y))")
            end
        end
    end

    println("动画录制完成！")
    println("输出文件: $output_path")

    # 清理资源
    destroy_scene(scene_ptr)
end

# 运行模拟
try
    simulation()
catch e
    println("发生错误: $e")
    println("错误堆栈:")
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
end