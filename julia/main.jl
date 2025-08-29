using Makie
using GLMakie
using Libdl
using Pkg
using FileIO

# 获取 Julia 脚本的目录
const SCRIPT_DIR = @__DIR__

# 构建一个绝对路径来定位 Rust 库文件
const LIB_PATH = abspath(joinpath(SCRIPT_DIR, "../target/release/libfluid_sim.so"))

if !isfile(LIB_PATH)
    error("找不到编译后的 Rust 库文件: $(LIB_PATH)。请确保你已成功运行 `cargo build --release`。")
end

const LIB = LIB_PATH

function setup_scene(num_particles::Int32)::Ptr{Cvoid}
    ccall((:setup_scene_wrapper, LIB), Ptr{Cvoid}, (Cint,), num_particles)
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

function get_output_dims()::Tuple{Int, Int}
    rows = ccall((:get_output_rows, LIB), Cint, ())
    cols = ccall((:get_output_cols, LIB), Cint, ())
    return Int(rows), Int(cols)
end

function get_output(scene_ptr::Ptr{Cvoid})
    rows, cols = get_output_dims()
    ptr = ccall((:get_output_wrapper, LIB), Ptr{Bool}, (Ptr{Cvoid},), scene_ptr)
    output_array = unsafe_wrap(Array, ptr, (rows, cols), own=true)
    return output_array
end

function main_simulation()
    
    NUM_PARTICLES = 400
    scene_ptr = setup_scene(Int32(NUM_PARTICLES))
    println("场景已初始化，粒子数量: $(NUM_PARTICLES)")
    
    set_gravity(scene_ptr, 0.0f0, -9.8f0)

    initial_data = get_output(scene_ptr)
    fluid_data = Observable(initial_data)

    set_theme!(backgroundcolor = :black)
    fig = Figure(resolution = (800, 800))
    ax = Axis(fig[1, 1], aspect = 1, yreversed = true)
    
    heatmap!(ax, fluid_data, colormap = [:black, :dodgerblue4], colorrange=(0,1))
    
    display(fig)

    for i in 1:200
        simulate_step(scene_ptr)
        new_data = get_output(scene_ptr)
        fluid_data[] = new_data
        
        sleep(0.01)
    end
    
    destroy_scene(scene_ptr)
end

main_simulation()