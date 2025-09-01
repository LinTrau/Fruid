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
    return Int(rows),
    Int(cols)
end

function get_output(scene_ptr::Ptr{Cvoid})
    rows, cols = get_output_dims(scene_ptr)
    ptr = ccall((:get_output_wrapper, LIB), Ptr{Bool}, (Ptr{Cvoid},), scene_ptr)
    output_array = unsafe_wrap(Array, ptr, (cols, rows), own=true)
    return output_array
end

function get_density(scene_ptr::Ptr{Cvoid})
    rows, cols = get_output_dims(scene_ptr)
    ptr = ccall((:get_density_wrapper, LIB), Ptr{Cfloat}, (Ptr{Cvoid},), scene_ptr)
    density_array = unsafe_wrap(Array, ptr, (cols, rows), own=false)
    return density_array
end

function get_u_velocity(scene_ptr::Ptr{Cvoid})
    rows, cols = get_output_dims(scene_ptr)
    ptr = ccall((:get_u_velocity_wrapper, LIB), Ptr{Cfloat}, (Ptr{Cvoid},), scene_ptr)
    u_array = unsafe_wrap(Array, ptr, (cols, rows), own=false)
    return u_array
end

function get_v_velocity(scene_ptr::Ptr{Cvoid})
    rows, cols = get_output_dims(scene_ptr)
    ptr = ccall((:get_v_velocity_wrapper, LIB), Ptr{Cfloat}, (Ptr{Cvoid},), scene_ptr)
    v_array = unsafe_wrap(Array, ptr, (cols, rows), own=false)
    return v_array
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

    gravity_x_str = replace(string(gravity[1]), "." => "p")
    gravity_y_str = replace(string(gravity[2]), "." => "p")
    filename = "particles_$(num_particles)_steps_$(sim_steps)_gravity_$(gravity_x_str)_$(gravity_y_str).mp4"
    output_path = joinpath(output_dir, filename)

    scene_ptr = setup_scene(Int32(num_particles), Int32(window_height), Int32(window_width))

    println("场景已初始化。")

    set_gravity(scene_ptr, Float32(gravity[1]), Float32(gravity[2]))

    # --- 更新 Makie 可视化设置 ---

    # 密度场数据
    fluid_data = Observable(get_density(scene_ptr))

    # 速度场数据
    u_vel = Observable(get_u_velocity(scene_ptr))
    v_vel = Observable(get_v_velocity(scene_ptr))

    # 粒子位置数据
    particles = Observable(get_particles(scene_ptr))
    particle_x = @lift($particles[1, :])
    particle_y = @lift($particles[2, :])

    set_theme!(backgroundcolor=:black)
    fig = Figure(resolution=(window_width, window_height))
    ax = Axis(fig[1, 1], aspect=DataAspect(), yreversed=true, title="Fluid simulation")

    # 1. 密度场 (背景)
    density_heatmap = heatmap!(ax, fluid_data, colormap=:dense)
    Colorbar(fig[1, 2], density_heatmap, label="Density")

    # 2. 速度场 (箭头)
    # 为了防止画面过于混乱，我们每隔10个点绘制一个箭头
    skip = 10
    rows, cols = get_output_dims(scene_ptr)
    x_pos = 1:skip:cols
    y_pos = 1:skip:rows
    # 使用 lift 保证箭头数据在每一帧都能更新
    u_points = @lift($u_vel[x_pos, y_pos])
    v_points = @lift($v_vel[x_pos, y_pos])
    arrows2d!(ax, x_pos, y_pos, u_points, v_points,
        tipwidth=7, lengthscale=0.3, color=:white)

    # 3. 粒子位置 (散点)
    # 你可以选择速度场或粒子位置，或者都显示。如果都显示，粒子图可能会覆盖速度场。
    # 这里我们默认显示粒子，如果你想看速度场，可以注释掉下面这行。
    scatter!(ax, particle_x, particle_y, color=:cyan, markersize=2)


    display(fig)

    record(fig, output_path, 1:sim_steps;
        framerate=30) do frame
        simulate_step(scene_ptr)

        # 在每一帧更新所有数据
        fluid_data[] = get_density(scene_ptr)
        u_vel[] = get_u_velocity(scene_ptr)
        v_vel[] = get_v_velocity(scene_ptr)
        particles[] = get_particles(scene_ptr)

        println("帧: ", frame)
    end
end

simulation()