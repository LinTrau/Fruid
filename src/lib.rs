#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(non_upper_case_globals)]

use libm::{floorf, sqrtf};
use std::ptr;

#[derive(PartialEq, Clone, Copy)]
pub enum CellType {
    AIR_CELL,
    FLUID_CELL,
    SOLID_CELL,
}

pub struct FlipFluid {
    density: f32,
    fNumX: f32,
    fNumY: f32,
    h: f32,
    fInvSpacing: f32,
    fNumCells: f32,
    pub u: Vec<f32>,
    pub v: Vec<f32>,
    du: Vec<f32>,
    dv: Vec<f32>,
    prevU: Vec<f32>,
    prevV: Vec<f32>,
    p: Vec<f32>,
    s: Vec<f32>,
    cellType: Vec<CellType>,
    pub particlePos: Vec<f32>,
    particleVel: Vec<f32>,
    particleDensity: Vec<f32>,
    particleRestDensity: f32,
    particleRadius: f32,
    pInvSpacing: f32,
    pNumX: i32,
    pNumY: i32,
    pNumCells: i32,
    numCellParticles: Vec<i32>,
    firstCellParticle: Vec<i32>,
    cellParticleIds: Vec<i32>,
    pub numParticles: i32,
    pub dens: Vec<f32>,
    pub dens_prev: Vec<f32>,
}

impl FlipFluid {
    fn new(
        density: f32,
        width: f32,
        height: f32,
        spacing: f32,
        particleRadius: f32,
        max_particles: i32,
        _num_rows: i32,
        _num_cols: i32,
    ) -> FlipFluid {
        let fNumX = floorf(width / spacing);
        let fNumY = floorf(height / spacing);
        let h = (width / fNumX).max(height / fNumY);
        let fInvSpacing = 1.0 / h;
        let fNumCells = fNumX * fNumY;
        let number_of_cells_setting = (fNumX * fNumY) as usize;
        let number_of_cells_x2_setting = number_of_cells_setting * 2;
        let number_of_cells_setting_plus1 = number_of_cells_setting + 1;
        let max_particles_usize = max_particles as usize;

        let u = vec![0.0; number_of_cells_x2_setting];
        let v = vec![0.0; number_of_cells_x2_setting];
        let du = vec![0.0; number_of_cells_x2_setting];
        let dv = vec![0.0; number_of_cells_x2_setting];
        let prevU = vec![0.0; number_of_cells_x2_setting];
        let prevV = vec![0.0; number_of_cells_x2_setting];
        let p = vec![0.0; number_of_cells_setting];
        let mut s = vec![0.0; number_of_cells_setting];
        let n = fNumY as usize;

        for i in 0..fNumX as usize {
            for j in 0..fNumY as usize {
                let mut s_val = 1.0; 
                if i == 0 || i == (fNumX as usize) - 1 || j == 0 || j == (fNumY as usize) - 1 {
                    s_val = 0.0;
                }
                s[i * n + j] = s_val;
            }
        }

        let cellType = vec![CellType::AIR_CELL; number_of_cells_setting];
        let mut particlePos = vec![0.0; max_particles_usize * 2];
        let mut count: usize = 0;

        let particle_fill_x_max = (max_particles as f32).sqrt().ceil() as i32;
        let particle_fill_y_max = (max_particles as f32).sqrt().ceil() as i32;

        for i in 1..=particle_fill_y_max {
            for j in 1..=particle_fill_x_max {
                if count < max_particles_usize {
                    particlePos[count * 2] = (j as f32) / 2.0;
                    particlePos[count * 2 + 1] = (i as f32) / 2.0;
                    count += 1;
                } else {
                    break;
                }
            }
        }
        let particleVel = vec![0.0; max_particles_usize * 2];
        let particleDensity = vec![0.0; number_of_cells_setting];
        let particleRestDensity = 0.0f32;

        let pInvSpacing = 1.0;
        let pNumX = floorf(width * pInvSpacing) as i32;
        let pNumY = floorf(height * pInvSpacing) as i32;
        let pNumCells = pNumX * pNumY;
        let numCellParticles = vec![0; number_of_cells_setting];
        let firstCellParticle = vec![0; number_of_cells_setting_plus1];
        let cellParticleIds = vec![0; max_particles_usize];
        let numParticles = 0i32;

        let dens = vec![0.0; number_of_cells_setting];
        let dens_prev = vec![0.0; number_of_cells_setting];

        FlipFluid {
            density,
            fNumX,
            fNumY,
            h,
            fInvSpacing,
            fNumCells,
            u,
            v,
            du,
            dv,
            prevU,
            prevV,
            p,
            s,
            cellType,
            particlePos,
            particleVel,
            particleDensity,
            particleRestDensity,
            particleRadius,
            pInvSpacing,
            pNumY,
            numCellParticles,
            cellParticleIds,
            pNumX,
            pNumCells,
            firstCellParticle,
            numParticles,
            dens,
            dens_prev,
        }
    }

    //--------------------------------------------------------------------------------
    // REFACTORED SIMULATION FLOW
    //--------------------------------------------------------------------------------
    
    /// FLIP模拟算法
    pub fn simulate(
        &mut self,
        dt: f32,
        x_gravity: f32,
        y_gravity: f32,
        flip_ratio: f32,
        num_pressure_iters: i32,
        num_particle_iters: i32,
        over_relaxation: f32,
        compensate_drift: bool,
        separate_particles: bool,
    ) {
        // 步骤 1: 将重力等外力施加到粒子速度上
        self.add_gravity_to_particles(dt, y_gravity, x_gravity);

        // 步骤 2: 通过推开粒子来防止粒子聚集，改善模拟效果
        if separate_particles {
            self.push_particles_apart(num_particle_iters);
        }

        // 步骤 3: 将粒子速度传递到网格 (PIC/FLIP的第一部分)。
        // 这个函数会保存压力求解前的速度场(prevU, prevV)。
        self.transfer_velocities(true, flip_ratio);

        // 步骤 4: 求解压力，强制网格速度场不可压缩。
        self.solve_incompressibility(num_pressure_iters, dt, over_relaxation, compensate_drift);

        // 步骤 5: 将网格速度的变化更新回粒子 (FLIP/PIC的第二部分)。
        self.transfer_velocities(false, flip_ratio);

        // 步骤 6: 使用更新后的粒子速度来移动粒子（平流）。
        self.advect_particles(dt);

        // 步骤 7: 处理粒子与容器边界的碰撞。
        self.handle_particle_collisions();

        // 步骤 8: 更新用于可视化的密度场。
        self.update_density();
    }
    
    //--------------------------------------------------------------------------------
    // CORE SIMULATION SUB-STEPS
    //--------------------------------------------------------------------------------

    /// 步骤 1: 施加重力到每个粒子上
    fn add_gravity_to_particles(&mut self, dt: f32, y_gravity: f32, x_gravity: f32) {
        for i in 0..self.numParticles as usize {
            self.particleVel[2 * i] += dt * x_gravity;
            self.particleVel[2 * i + 1] += dt * y_gravity;
        }
    }

    /// 步骤 2: 将过于靠近的粒子推开，防止它们聚集在一起。
    fn push_particles_apart(&mut self, numIters: i32) {
        // ... (此函数实现完整，无需修改)
        self.numCellParticles.fill(0);

        for i in 0..self.numParticles as usize {
            let x: f32 = self.particlePos[2 * i];
            let y: f32 = self.particlePos[2 * i + 1];

            let xi = clamp(floorf(x) as i32, 1, self.pNumX - 2);
            let yi = clamp(floorf(y) as i32, 1, self.pNumY - 2);
            let celNr = xi * self.pNumY + yi;

            self.numCellParticles[celNr as usize] += 1;
        }

        let mut first = 0;
        for i in 0..self.pNumCells as usize {
            first += self.numCellParticles[i];
            self.firstCellParticle[i] = first;
        }
        self.firstCellParticle[self.pNumCells as usize] = first;

        for i in 0..self.numParticles as usize {
            let x = self.particlePos[2 * i];
            let y = self.particlePos[2 * i + 1];

            let xi = clamp(floorf(x * self.pInvSpacing) as i32, 1, self.pNumX - 2);
            let yi = clamp(floorf(y * self.pInvSpacing) as i32, 1, self.pNumY - 2);
            let cellNr = xi * self.pNumY + yi;

            self.firstCellParticle[cellNr as usize] -= 1;
            self.cellParticleIds[self.firstCellParticle[cellNr as usize] as usize] = i as i32;
        }

        let minDist = 2.0 * self.particleRadius;
        let minDist2 = minDist * minDist;

        for _i in 0..numIters {
            for i in 0..self.numParticles as usize {
                let px = self.particlePos[2 * i];
                let py = self.particlePos[2 * i + 1];

                let pxi = floorf(px * self.pInvSpacing) as i32;
                let pyi = floorf(py * self.pInvSpacing) as i32;
                let x0 = (pxi - 1).max(0);
                let y0 = (pyi - 1).max(0);
                let x1 = (pxi + 1).min(self.pNumX - 1);
                let y1 = (pyi + 1).min(self.pNumY - 1);

                for xi in x0..=x1 {
                    for yi in y0..=y1 {
                        let cellNr = xi * self.pNumY + yi;
                        let first = self.firstCellParticle[cellNr as usize];
                        let last = self.firstCellParticle[(cellNr + 1) as usize];
                        for j in first..last {
                            let id = self.cellParticleIds[j as usize];
                            if id as usize == i {
                                continue;
                            }

                            let qx = self.particlePos[(2 * id) as usize];
                            let qy = self.particlePos[(2 * id + 1) as usize];

                            let mut dx = qx - px;
                            let mut dy = qy - py;
                            let d2 = dx * dx + dy * dy;

                            if d2 > minDist2 || d2 == 0f32 {
                                continue;
                            }
                            let d = sqrtf(d2);
                            let s = 0.5f32 * (minDist - d) / d;
                            dx *= s;
                            dy *= s;
                            self.particlePos[2 * i] -= dx;
                            self.particlePos[2 * i + 1] -= dy;
                            self.particlePos[(2 * id) as usize] += dx;
                            self.particlePos[(2 * id + 1) as usize] += dy;
                        }
                    }
                }
            }
        }
    }
    
    /// 步骤 3 & 5: 在粒子和网格之间传递速度。
    fn transfer_velocities(&mut self, toGrid: bool, flipRatio: f32) {
        let n = self.fNumY;
        let h = self.h;
        let h1 = self.fInvSpacing;
        let h2 = 0.5 * h;

        if toGrid {
            self.prevU = self.u.clone();
            self.prevV = self.v.clone();
            self.du.fill(0.0);
            self.dv.fill(0.0);
            self.u.fill(0.0);
            self.v.fill(0.0);

            for i in 0..self.fNumCells as usize {
                self.cellType[i] = if self.s[i] == 0.0 {
                    CellType::SOLID_CELL
                } else {
                    CellType::AIR_CELL
                };
            }

            for i in 0..self.numParticles as usize {
                let x = self.particlePos[2 * i];
                let y = self.particlePos[2 * i + 1];

                let xi = clamp(floorf(x * h1), 1.0, self.fNumX - 2.0) as usize;
                let yi = clamp(floorf(y * h1), 1.0, self.fNumY - 2.0) as usize;
                let cellNr = xi * n as usize + yi;

                if self.cellType[cellNr] == CellType::AIR_CELL {
                    self.cellType[cellNr] = CellType::FLUID_CELL;
                }
            }
        }

        for component in 0..2 {
            let dx = if component == 0 { 0.0 } else { h2 };
            let dy = if component == 0 { h2 } else { 0.0 };
            let f = if component == 0 {
                &mut self.u
            } else {
                &mut self.v
            };
            let prevF = if component == 0 {
                &self.prevU
            } else {
                &self.prevV
            };
            let d = if component == 0 {
                &mut self.du
            } else {
                &mut self.dv
            };

            for i in 0..self.numParticles as usize {
                let mut x = self.particlePos[2 * i];
                let mut y = self.particlePos[2 * i + 1];

                x = clamp(x, h, (self.fNumX - 1.0) * h);
                y = clamp(y, h, (self.fNumY - 1.0) * h);

                let x0 = floorf((x - dx) * h1).min(self.fNumX - 2.0) as usize;
                let tx = ((x - dx) - x0 as f32 * h) * h1;
                let x1 = ((x0 as f32) + 1.0).min(self.fNumX - 2.0) as usize;

                let y0 = floorf((y - dy) * h1).min(self.fNumY - 2.0) as usize;
                let ty = ((y - dy) - y0 as f32 * h) * h1;
                let y1 = ((y0 as f32) + 1.0).min(self.fNumY - 2.0) as usize;

                let sx = 1.0 - tx;
                let sy = 1.0 - ty;
                let d0 = sx * sy;
                let d1 = tx * sy;
                let d2 = tx * ty;
                let d3 = sx * ty;

                let nr0 = x0 * n as usize + y0;
                let nr1 = x1 * n as usize + y0;
                let nr2 = x1 * n as usize + y1;
                let nr3 = x0 * n as usize + y1;

                if toGrid {
                    let pv = self.particleVel[2 * i + component];
                    f[nr0] += pv * d0;
                    d[nr0] += d0;
                    f[nr1] += pv * d1;
                    d[nr1] += d1;
                    f[nr2] += pv * d2;
                    d[nr2] += d2;
                    f[nr3] += pv * d3;
                    d[nr3] += d3;
                } else {
                    let offset = if component == 0 { n } else { 1.0 };
                    let valid0 = if self.cellType[nr0] != CellType::AIR_CELL
                        || self.cellType[nr0 - offset as usize] != CellType::AIR_CELL
                    {
                        1.0
                    } else {
                        0.0
                    };
                    let valid1 = if self.cellType[nr1] != CellType::AIR_CELL
                        || self.cellType[nr1 - offset as usize] != CellType::AIR_CELL
                    {
                        1.0
                    } else {
                        0.0
                    };
                    let valid2 = if self.cellType[nr2] != CellType::AIR_CELL
                        || self.cellType[nr2 - offset as usize] != CellType::AIR_CELL
                    {
                        1.0
                    } else {
                        0.0
                    };
                    let valid3 = if self.cellType[nr3] != CellType::AIR_CELL
                        || self.cellType[nr3 - offset as usize] != CellType::AIR_CELL
                    {
                        1.0
                    } else {
                        0.0
                    };

                    let v = self.particleVel[2 * i + component];
                    let d_val = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

                    if d_val > 0.0 {
                        let picV = (valid0 * d0 * f[nr0]
                            + valid1 * d1 * f[nr1]
                            + valid2 * d2 * f[nr2]
                            + valid3 * d3 * f[nr3])
                            / d_val;
                        let corr = (valid0 * d0 * (f[nr0] - prevF[nr0])
                            + valid1 * d1 * (f[nr1] - prevF[nr1])
                            + valid2 * d2 * (f[nr2] - prevF[nr2])
                            + valid3 * d3 * (f[nr3] - prevF[nr3]))
                            / d_val;
                        let flipV = v + corr;

                        self.particleVel[2 * i + component] =
                            (1.0 - flipRatio) * picV + flipRatio * flipV;
                    }
                }
            }

            if toGrid {
                for i in 0..f.len() {
                    if d[i] > 0.0 {
                        f[i] /= d[i];
                    }
                }

                for i in 0..self.fNumX as usize {
                    for j in 0..self.fNumY as usize {
                        let solid = self.cellType[i * n as usize + j] == CellType::SOLID_CELL;
                        if solid
                            || (i > 0
                                && self.cellType[(i - 1) * n as usize + j] == CellType::SOLID_CELL)
                        {
                            self.u[i * n as usize + j] = self.prevU[i * n as usize + j];
                        }
                        if solid
                            || (j > 0
                                && self.cellType[i * n as usize + j - 1] == CellType::SOLID_CELL)
                        {
                            self.v[i * n as usize + j] = self.prevV[i * n as usize + j];
                        }
                    }
                }
            }
        }
    }
    
    /// 步骤 4: 这是一个完整的压力求解和应用函数，确保流体不可压缩。
    fn solve_incompressibility(
        &mut self,
        numIters: i32,
        dt: f32,
        overRelaxation: f32,
        compensateDrift: bool,
    ) {
        self.p.fill(0.0);
        self.prevU = self.u.clone();
        self.prevV = self.v.clone();

        let n = self.fNumY;
        let cp = self.density * self.h / dt;

        for _iter in 0..numIters {
            for i in 1..self.fNumX as usize {
                for j in 1..self.fNumY as usize {
                    if self.cellType[i * n as usize + j] != CellType::FLUID_CELL {
                        continue;
                    }

                    let center = i * n as usize + j;
                    let left = (i - 1) * n as usize + j;
                    let right = (i + 1) * n as usize + j;
                    let bottom = i * n as usize + j - 1;
                    let top = i * n as usize + j + 1;

                    let _s = self.s[center];
                    let sx0 = self.s[left];
                    let sx1 = self.s[right];
                    let sy0 = self.s[bottom];
                    let sy1 = self.s[top];
                    let s_val = sx0 + sx1 + sy0 + sy1;
                    if s_val == 0.0 {
                        continue;
                    }

                    let mut div = self.u[right] - self.u[center] + self.v[top] - self.v[center];

                    if self.particleRestDensity > 0.0 && compensateDrift {
                        let k = 1.0;
                        let compression =
                            self.particleDensity[i * n as usize + j] - self.particleRestDensity;
                        if compression > 0.0 {
                            div -= k * compression;
                        }
                    }

                    let mut p = -div / s_val;
                    p *= overRelaxation;
                    self.p[center] += cp * p;

                    self.u[center] -= sx0 * p;
                    self.u[right] += sx1 * p;
                    self.v[center] -= sy0 * p;
                    self.v[top] += sy1 * p;
                }
            }
        }
    }

    /// 步骤 6: 使用最终速度移动粒子
    fn advect_particles(&mut self, dt: f32) {
        for i in 0..self.numParticles as usize {
            self.particlePos[i * 2] += self.particleVel[i * 2] * dt;
            self.particlePos[i * 2 + 1] += self.particleVel[i * 2 + 1] * dt;
        }
    }
    
    /// 步骤 7: 粒子与边界的碰撞处理
    fn handle_particle_collisions(&mut self) {
        let minX = 1.0;
        let maxX = self.fNumX - 1.0;
        let minY = 1.0;
        let maxY = self.fNumY - 1.0;

        for i in 0..self.numParticles as usize {
            let mut x = self.particlePos[2 * i];
            let mut y = self.particlePos[2 * i + 1];

            if x < minX {
                x = minX;
                self.particleVel[2 * i] = 0.0;
            }
            if x > maxX {
                x = maxX;
                self.particleVel[2 * i] = 0.0;
            }
            if y < minY {
                y = minY;
                self.particleVel[2 * i + 1] = 0.0;
            }
            if y > maxY {
                y = maxY;
                self.particleVel[2 * i + 1] = 0.0;
            }
            self.particlePos[2 * i] = x;
            self.particlePos[2 * i + 1] = y;
        }
    }
    
    /// 步骤 8: 从粒子位置更新密度场
    pub fn update_density(&mut self) {
        self.dens.fill(0.0);
        for i in 0..self.numParticles as usize {
            let x = self.particlePos[i * 2];
            let y = self.particlePos[i * 2 + 1];

            if x < 0.0 || x >= self.fNumX || y < 0.0 || y >= self.fNumY { continue; }

            let i_cell = floorf(x) as usize;
            let j_cell = floorf(y) as usize;

            if i_cell < self.fNumX as usize -1 && j_cell < self.fNumY as usize -1 {
                let x_frac = x - i_cell as f32;
                let y_frac = y - j_cell as f32;

                let w00 = (1.0 - x_frac) * (1.0 - y_frac);
                let w10 = x_frac * (1.0 - y_frac);
                let w01 = (1.0 - x_frac) * y_frac;
                let w11 = x_frac * y_frac;
                
                let idx00 = j_cell * self.fNumX as usize + i_cell;
                let idx10 = j_cell * self.fNumX as usize + i_cell + 1;
                let idx01 = (j_cell + 1) * self.fNumX as usize + i_cell;
                let idx11 = (j_cell + 1) * self.fNumX as usize + i_cell + 1;

                self.dens[idx00] += w00;
                self.dens[idx10] += w10;
                self.dens[idx01] += w01;
                self.dens[idx11] += w11;
            }
        }
    }

    //--------------------------------------------------------------------------------
    // UTILITY AND OTHER FUNCTIONS 
    //--------------------------------------------------------------------------------

    pub fn add_initial_fluid(&mut self) {
        let start_x = self.fNumX / 4.0;
        let end_x = self.fNumX / 2.0;
        let start_y = self.fNumY / 4.0;
        let end_y = self.fNumY / 2.0;

        let mut current_particle_count = 0;
        for x in (start_x as usize)..(end_x as usize) {
            for y in (start_y as usize)..(end_y as usize) {
                if current_particle_count < self.numParticles as usize {
                    self.particlePos[current_particle_count * 2] = x as f32;
                    self.particlePos[current_particle_count * 2 + 1] = y as f32;
                    current_particle_count += 1;
                } else {
                    break;
                }
            }
            if current_particle_count >= self.numParticles as usize {
                break;
            }
        }

        self.numParticles = current_particle_count as i32;
    }

    /// [INFO] 下面这些函数用于经典的欧拉烟雾/染料模拟。
    /// 它们在当前的FLIP模拟流程中没有被调用，但如果你想在流体中模拟烟雾的扩散和流动，可以在 `simulate` 函数的末尾调用 `dens_step`。
    
    pub fn add_density_source(&mut self, source_x: i32, source_y: i32, amount: f32) {
        if source_x >= 0
            && source_x < self.fNumX as i32
            && source_y >= 0
            && source_y < self.fNumY as i32
        {
            self.dens_prev[source_x as usize * self.fNumY as usize + source_y as usize] += amount;
        }
    }

    fn dens_step(&mut self, dt: f32, diff: f32) {
        let mut x = self.dens.clone();
        let mut x0 = self.dens_prev.clone();
        let u = self.u.clone();
        let v = self.v.clone();
        self.diffuse(1, &mut x0, &mut x, diff, dt, self.fNumY as i32);
        self.advect(
            1,
            &mut x,
            &x0,
            &u,
            &v,
            dt,
            self.fNumX as i32,
            self.fNumY as i32,
        );
        self.dens = x;
        self.dens_prev.fill(0.0);
    }
    
    fn diffuse(&self, _b: i32, x: &mut Vec<f32>, x0: &Vec<f32>, diff: f32, dt: f32, N: i32) {
        let a = dt * diff * (N - 2) as f32 * (N - 2) as f32;
        for _k in 0..20 {
            for i in 1..N - 1 {
                for j in 1..N - 1 {
                    x[(i * N + j) as usize] = (x0[(i * N + j) as usize]
                        + a * (x[((i + 1) * N + j) as usize]
                            + x[((i - 1) * N + j) as usize]
                            + x[(i * N + j + 1) as usize]
                            + x[(i * N + j - 1) as usize]))
                        / (1.0 + 4.0 * a);
                }
            }
            self.set_bnd(_b, x, N);
        }
    }

    fn advect(
        &self,
        _b: i32,
        d: &mut Vec<f32>,
        d0: &Vec<f32>,
        u: &Vec<f32>,
        v: &Vec<f32>,
        dt: f32,
        N: i32,
        N_y: i32,
    ) {
        let N_f32 = N as f32;
        let dt0 = dt * N_f32;
        for i in 1..N - 1 {
            for j in 1..N_y - 1 {
                let mut x = (i as f32) - dt0 * u[(i * N_y + j) as usize];
                let mut y = (j as f32) - dt0 * v[(i * N_y + j) as usize];
                if x < 0.5 { x = 0.5; }
                if x > N_f32 + 0.5 { x = N_f32 + 0.5; }
                let i0 = floorf(x);
                let i1 = i0 + 1.0;
                if y < 0.5 { y = 0.5; }
                if y > N_y as f32 + 0.5 { y = N_y as f32 + 0.5; }
                let j0 = floorf(y);
                let j1 = j0 + 1.0;
                let s1 = x - i0;
                let s0 = 1.0 - s1;
                let t1 = y - j0;
                let t0 = 1.0 - t1;
                d[(i * N_y + j) as usize] = s0
                    * (t0 * d0[i0 as usize * N_y as usize + j0 as usize]
                        + t1 * d0[i0 as usize * N_y as usize + j1 as usize])
                    + s1 * (t0 * d0[i1 as usize * N_y as usize + j0 as usize]
                        + t1 * d0[i1 as usize * N_y as usize + j1 as usize]);
            }
        }
        self.set_bnd(_b, d, N_y);
    }
    
    fn set_bnd(&self, b: i32, x: &mut Vec<f32>, N: i32) {
        for i in 1..N - 1 {
            x[(i * N) as usize] = if b == 2 { -x[(i * N + 1) as usize] } else { x[(i * N + 1) as usize] };
            x[(i * N + N - 1) as usize] = if b == 2 { -x[(i * N + N - 2) as usize] } else { x[(i * N + N - 2) as usize] };
        }
        for j in 1..N - 1 {
            x[j as usize] = if b == 1 { -x[(1 * N + j) as usize] } else { x[(1 * N + j) as usize] };
            x[((N - 1) * N + j) as usize] = if b == 1 { -x[((N - 2) * N + j) as usize] } else { x[((N - 2) * N + j) as usize] };
        }
        x[0] = 0.5 * (x[1] + x[N as usize]);
        x[(N - 1) as usize] = 0.5 * (x[(N - 2) as usize] + x[(2 * N - 1) as usize]);
        x[((N - 1) * N) as usize] = 0.5 * (x[((N - 1) * N + 1) as usize] + x[((N - 2) * N) as usize]);
        x[((N - 1) * N + N - 1) as usize] = 0.5 * (x[((N - 1) * N + N - 2) as usize] + x[((N - 2) * N + N - 1) as usize]);
    }
}

pub struct Scene {
    xGravity: f32,
    yGravity: f32,
    dt: f32,
    flipRatio: f32,
    numPressureIters: i32,
    numParticleIters: i32,
    frameNr: i32,
    overRelaxation: f32,
    compensateDrift: bool,
    separateParticles: bool,
    paused: bool,
    pub fluid: FlipFluid,
    num_cols: i32,
    num_rows: i32,
}

impl Scene {
    pub fn setup_scene(num_particles: i32, num_rows: i32, num_cols: i32) -> Scene {
        let xGravity = 0.0;
        let yGravity = -9.81; 
        let flipRatio = 0.9;
        let frameNr = 0;
        let overRelaxation = 1.9;
        let compensateDrift = true;
        let separateParticles = false;
        let paused = false;
        let dt = 1.0 / 60.0;
        let numPressureIters = 50; 
        let numParticleIters = 2;

        let res_width = num_cols as f32;
        let res_height = num_rows as f32;

        let tankHeight = 1.0 * res_height;
        let tankWidth = 1.0 * res_width;
        let h = tankHeight / res_height;
        let density = 1000.0;
        let r = 0.3 * h; 

        let mut fluid = FlipFluid::new(
            density,
            tankWidth,
            tankHeight,
            h,
            r,
            num_particles,
            num_rows,
            num_cols,
        );
        fluid.numParticles = num_particles;

        Scene {
            xGravity,
            yGravity,
            dt,
            flipRatio,
            numPressureIters,
            numParticleIters,
            frameNr,
            overRelaxation,
            compensateDrift,
            separateParticles,
            paused,
            fluid,
            num_cols,
            num_rows,
        }
    }
    
    pub fn pause(&mut self) { self.paused = true; }
    pub fn unpause(&mut self) { self.paused = false; }
    pub fn is_paused(&self) -> bool { self.paused }

    pub fn simulate(&mut self) {
        if self.paused { return; }
        
        self.fluid.simulate(
            self.dt,
            self.xGravity,
            self.yGravity,
            self.flipRatio,
            self.numPressureIters,
            self.numParticleIters,
            self.overRelaxation,
            self.compensateDrift,
            self.separateParticles,
        );
        self.frameNr += 1;
    }
    
    pub fn set_gravity(&mut self, accel_measurment: [f32; 2]) {
        self.xGravity = accel_measurment[0];
        self.yGravity = accel_measurment[1];
    }
    
    pub fn get_output(&mut self) -> Vec<bool> {
        let mut output_frame: Vec<bool> = vec![false; (self.num_rows * self.num_cols) as usize];
        let n = self.fluid.fNumY as usize;
        for i in 0..self.num_rows as usize {
            for j in 0..self.num_cols as usize {
                if j < self.fluid.fNumX as usize && i < self.fluid.fNumY as usize {
                    if self.fluid.cellType[j * n + i] == CellType::FLUID_CELL {
                        output_frame[i * self.num_cols as usize + j] = true;
                    }
                }
            }
        }
        output_frame
    }

    pub fn get_density_data(&mut self) -> *const f32 {
        self.fluid.dens.as_ptr()
    }
}

//--------------------------------------------------------------------------------
// Julia FFI WRAPPERS 
//--------------------------------------------------------------------------------

#[unsafe(no_mangle)]
pub extern "C" fn setup_scene_wrapper(
    num_particles: i32,
    num_rows: i32,
    num_cols: i32,
) -> *mut Scene {
    let mut scene = Box::new(Scene::setup_scene(num_particles, num_rows, num_cols));
    scene.fluid.add_initial_fluid();
    Box::into_raw(scene)
}

#[unsafe(no_mangle)]
pub extern "C" fn destroy_scene_wrapper(ptr: *mut Scene) {
    if ptr.is_null() { return; }
    unsafe { drop(Box::from_raw(ptr)); }
}

#[unsafe(no_mangle)]
pub extern "C" fn simulate_wrapper(ptr: *mut Scene) {
    if ptr.is_null() { return; }
    let scene = unsafe { &mut *ptr };
    scene.simulate();
}

#[unsafe(no_mangle)]
pub extern "C" fn set_gravity_wrapper(ptr: *mut Scene, x_gravity: f32, y_gravity: f32) {
    if ptr.is_null() { return; }
    let scene = unsafe { &mut *ptr };
    scene.set_gravity([x_gravity, y_gravity]);
}

#[unsafe(no_mangle)]
pub extern "C" fn get_output_wrapper(ptr: *mut Scene) -> *const bool {
    if ptr.is_null() { return ptr::null(); }
    let scene = unsafe { &mut *ptr };
    let flat_output = scene.get_output();
    let boxed_output = flat_output.into_boxed_slice();
    Box::into_raw(boxed_output) as *const bool
}

#[unsafe(no_mangle)]
pub extern "C" fn get_density_wrapper(ptr: *mut Scene) -> *const f32 {
    if ptr.is_null() { return ptr::null(); }
    let scene = unsafe { &mut *ptr };
    scene.get_density_data()
}

#[unsafe(no_mangle)]
pub extern "C" fn get_output_rows(ptr: *mut Scene) -> i32 {
    if ptr.is_null() { return 0; }
    let scene = unsafe { &*ptr };
    scene.num_rows
}

#[unsafe(no_mangle)]
pub extern "C" fn get_output_cols(ptr: *mut Scene) -> i32 {
    if ptr.is_null() { return 0; }
    let scene = unsafe { &*ptr };
    scene.num_cols
}

fn clamp<T: PartialOrd>(x: T, min: T, max: T) -> T {
    if x < min { min } else if x > max { max } else { x }
}