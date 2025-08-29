using Pkg

function install_dependencies()
    packages = ["TOML", "Makie", "GLMakie", "Libdl", "FileIO"]
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

const SCRIPT_DIR = @__DIR__
const LIB_PATH = abspath(joinpath(SCRIPT_DIR, "../target/release/libfluid_sim.so"))
const CONFIG_PATH = abspath(joinpath(SCRIPT_DIR, "../config.toml"))
const LIB = LIB_PATH

if !isfile(LIB_PATH)
    error("找不到编译后的 Rust 库文件: $(LIB_PATH)。请确保你已成功运行 `cargo build --release`。")
end

function read_config(config_file::String)
    return TOML.parsefile(config_file)
end

# Rust函数调用
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

function get_output(scene_ptr::Ptr{Cvoid}, rows::Int32, cols::Int32)
    ptr = ccall((:get_output_wrapper, LIB), Ptr{Bool}, (Ptr{Cvoid},), scene_ptr)
    output_array = unsafe_wrap(Array, ptr, (rows, cols), own=true)
    return output_array
end

function simulation(config_file::String)
    if !isfile(config_file)
        error("配置文件未找到: $config_file")
    end

    config = read_config(config_file)
    
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
    
    initial_data = get_output(scene_ptr, Int32(window_height), Int32(window_width))
    fluid_data = Observable(initial_data)

    set_theme!(backgroundcolor = :black)
    fig = Figure(resolution = (window_width, window_height))
    ax = Axis(fig[1, 1], aspect = DataAspect(), yreversed = true)
    
    heatmap!(ax, fluid_data, colormap = [:black, :dodgerblue4], colorrange=(0,1))
    
    display(fig)

    record(fig, output_path, 1:sim_steps) do i
        simulate_step(scene_ptr)
        new_data = get_output(scene_ptr, Int32(window_height), Int32(window_width))
        fluid_data[] = new_data
        
        sleep(0.01)
    end
    
    destroy_scene(scene_ptr)
end

if isempty(ARGS)
    simulation(CONFIG_PATH)
else
    simulation(ARGS[1])
end