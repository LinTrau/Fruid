#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(non_upper_case_globals)]

pub mod FluidSimulation {

    use libm::{floorf, sqrtf};
    use std::ptr;

    static max_particles_setting: usize = 800;
    static number_of_vertical_cells_setting: usize = 23;
    static number_of_horizontal_cells_setting: usize = 23;
    static max_particles_x2_setting: usize = max_particles_setting * 2;
    static number_of_cells_setting: usize =
        number_of_vertical_cells_setting * number_of_horizontal_cells_setting;
    static number_of_cells_x2_setting: usize = number_of_cells_setting * 2;
    static number_of_cells_setting_plus1: usize = number_of_cells_setting + 1;

    static simHeight: f32 = 23.0;
    static simWidth: f32 = 23.0;

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
        u: [f32; number_of_cells_x2_setting],
        v: [f32; number_of_cells_x2_setting],
        du: [f32; number_of_cells_x2_setting],
        dv: [f32; number_of_cells_x2_setting],
        prevU: [f32; number_of_cells_x2_setting],
        prevV: [f32; number_of_cells_x2_setting],
        p: [f32; number_of_cells_setting],
        s: [f32; number_of_cells_setting],
        cellType: [CellType; number_of_cells_setting],
        _maxParticles: i32,
        pub particlePos: [f32; max_particles_x2_setting],
        particleVel: [f32; max_particles_x2_setting],
        particleDensity: [f32; number_of_cells_setting],
        particleRestDensity: f32,
        particleRadius: f32,
        pInvSpacing: f32,
        pNumX: i32,
        pNumY: i32,
        pNumCells: i32,
        numCellParticles: [i32; number_of_cells_setting],
        firstCellParticle: [i32; number_of_cells_setting_plus1],
        cellParticleIds: [i32; max_particles_setting],
        pub numParticles: i32,
    }

    impl FlipFluid {
        fn new(
            density: f32,
            width: f32,
            height: f32,
            spacing: f32,
            particleRadius: f32,
            _maxParticles: i32,
        ) -> FlipFluid {
            let fNumX = floorf(width / spacing);
            let fNumY = floorf(height / spacing);
            let h = (width / fNumX).max(height / fNumY);
            let fInvSpacing = 1.0 / h;
            let fNumCells = fNumX * fNumY;

            let u = [0f32; number_of_cells_x2_setting];
            let v = [0f32; number_of_cells_x2_setting];
            let du = [0f32; number_of_cells_x2_setting];
            let dv = [0f32; number_of_cells_x2_setting];
            let prevU = [0f32; number_of_cells_x2_setting];
            let prevV = [0f32; number_of_cells_x2_setting];
            let p = [0f32; number_of_cells_setting];
            let s = [0f32; number_of_cells_setting];
            let cellType = [CellType::AIR_CELL; number_of_cells_setting];
            let mut particlePos = [0.0; max_particles_x2_setting];
            let mut count: usize = 0;
            for i in 1..21 {
                for j in 1..21 {
                    particlePos[count * 2] = (j as f32) / 2.0;
                    particlePos[count * 2 + 1] = (i as f32) / 2.0;
                    count += 1;
                }
            }
            let particleVel = [0.0; max_particles_x2_setting];
            let particleDensity = [0f32; number_of_cells_setting];
            let particleRestDensity = 0.0f32;

            let pInvSpacing = 1.0;
            let pNumX = floorf(width * pInvSpacing) as i32;
            let pNumY = floorf(height * pInvSpacing) as i32;
            let pNumCells = pNumX * pNumY;
            let numCellParticles = [0; number_of_cells_setting];
            let firstCellParticle = [0; number_of_cells_setting_plus1];
            let cellParticleIds = [0; max_particles_setting];
            let numParticles = 0i32;

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
                _maxParticles,
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
            }
        }

        fn integrateParticles(&mut self, dt: f32, yGravity: f32, xGravity: f32) {
            for i in 0..self.numParticles as usize {
                self.particleVel[2 * i] += dt * xGravity;
                self.particleVel[2 * i + 1] += dt * yGravity;
                self.particlePos[2 * i] += self.particleVel[2 * i] * dt;
                self.particlePos[2 * i + 1] += self.particleVel[2 * i + 1] * dt;
            }
        }

        fn showParticles(&mut self) {
            for i in 0..self.numParticles as usize {
                let cell_location_y = floorf(self.particlePos[2 * i]) as usize;
                let cell_location_x = floorf(self.particlePos[2 * i + 1]) as usize;

                if cell_location_y < 23 && cell_location_x < 23 {
                    self.cellType[cell_location_x * 23 + cell_location_y] = CellType::FLUID_CELL;
                }
            }
        }

        fn pushParticlesApart(&mut self, numIters: i32) {
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

        fn handleParticleCollisions(&mut self) {
            let minX = 1.0;
            let maxX = 21.0;
            let minY = 1.0;
            let maxY = 21.0;

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

        fn updateParticleDensity(&mut self) {
            let n = self.fNumY;
            let h = self.h;
            let h1 = self.fInvSpacing;
            let h2 = 0.5 * h;

            let mut d = self.particleDensity.clone();
            d.fill(0.0);

            for i in 0..self.numParticles as usize {
                let mut x = self.particlePos[2 * i];
                let mut y = self.particlePos[2 * i + 1];

                x = clamp(x, h, (self.fNumX - 1.0) * h);
                y = clamp(y, h, (self.fNumY - 1.0) * h);

                let x0 = floorf((x - h2) * h1) as usize;
                let tx = ((x - h2) - x0 as f32 * h) * h1;
                let x1 = ((x0 as f32) + 1.0).min(self.fNumX - 2.0) as usize;

                let y0 = floorf((y - h2) * h1) as usize;
                let ty = ((y - h2) - y0 as f32 * h) * h1;
                let y1 = ((y0 as f32) + 1.0).min(self.fNumY - 2.0) as usize;

                let sx = 1.0 - tx;
                let sy = 1.0 - ty;

                if (x0 as f32) < self.fNumX && (y0 as f32) < self.fNumY {
                    d[x0 * n as usize + y0] += sx * sy;
                }
                if (x1 as f32) < self.fNumX && (y0 as f32) < self.fNumY {
                    d[x1 * n as usize + y0] += tx * sy;
                }
                if (x1 as f32) < self.fNumX && (y1 as f32) < self.fNumY {
                    d[x1 * n as usize + y1] += tx * ty;
                }
                if (x0 as f32) < self.fNumX && (y1 as f32) < self.fNumY {
                    d[x0 * n as usize + y1] += sx * ty;
                }
            }

            if self.particleRestDensity == 0.0 {
                let mut sum = 0.0;
                let mut numFluidCells = 0;

                for i in 0..self.fNumCells as usize {
                    if self.cellType[i] == CellType::FLUID_CELL {
                        sum += d[i];
                        numFluidCells += 1;
                    }
                }

                if numFluidCells > 0 {
                    self.particleRestDensity = sum / numFluidCells as f32;
                }
            }
        }

        fn transferVelocities(&mut self, toGrid: bool, flipRatio: f32) {
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
                let mut f = if component == 0 {
                    self.u.clone()
                } else {
                    self.v.clone()
                };
                let prevF = if component == 0 {
                    self.prevU.clone()
                } else {
                    self.prevV.clone()
                };
                let mut d = if component == 0 {
                    self.du.clone()
                } else {
                    self.dv.clone()
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
                        let d = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

                        if d > 0.0 {
                            let picV = (valid0 * d0 * f[nr0]
                                + valid1 * d1 * f[nr1]
                                + valid2 * d2 * f[nr2]
                                + valid3 * d3 * f[nr3])
                                / d;
                            let corr = (valid0 * d0 * (f[nr0] - prevF[nr0])
                                + valid1 * d1 * (f[nr1] - prevF[nr1])
                                + valid2 * d2 * (f[nr2] - prevF[nr2])
                                + valid3 * d3 * (f[nr3] - prevF[nr3]))
                                / d;
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
                                    && self.cellType[(i - 1) * n as usize + j]
                                        == CellType::SOLID_CELL)
                            {
                                self.u[i * n as usize + j] = self.prevU[i * n as usize + j];
                            }
                            if solid
                                || (j > 0
                                    && self.cellType[i * n as usize + j - 1]
                                        == CellType::SOLID_CELL)
                            {
                                self.v[i * n as usize + j] = self.prevV[i * n as usize + j];
                            }
                        }
                    }
                }
            }
        }

        fn solveIncompressibility(
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

            for _i in 0..self.fNumCells as usize {}

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
                        let s = sx0 + sx1 + sy0 + sy1;
                        if s == 0.0 {
                            continue;
                        }

                        let mut div = self.u[right] - self.u[center] + self.v[top] - self.v[center];

                        if self.particleRestDensity > 0.0 && compensateDrift {
                            let k = 1.0;
                            let compression =
                                self.particleDensity[i * n as usize + j] - self.particleRestDensity;
                            if compression > 0.0 {
                                div = div - k * compression;
                            }
                        }

                        let mut p = -div / s;
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

        pub fn simulate(
            &mut self,
            dt: f32,
            xGravity: f32,
            yGravity: f32,
            flipRatio: f32,
            numPressureIters: i32,
            numParticleIters: i32,
            overRelaxation: f32,
            compensateDrift: bool,
            separateParticles: bool,
        ) {
            let numSubSteps = 1.0;
            let sdt = dt / numSubSteps;

            for _i in 0..numSubSteps as usize {
                self.integrateParticles(sdt, yGravity, xGravity);
                if separateParticles {
                    self.pushParticlesApart(numParticleIters);
                }
                self.handleParticleCollisions();
                self.pushParticlesApart(numParticleIters);
                self.handleParticleCollisions();
                self.transferVelocities(true, 1.9);
                self.updateParticleDensity();
                self.solveIncompressibility(numPressureIters, sdt, overRelaxation, compensateDrift);
                self.transferVelocities(false, flipRatio);
                self.showParticles();
            }
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
    }

    impl Scene {
        pub fn setupScene(particles: i32) -> Scene {
            let xGravity = 0.0;
            let yGravity = 0.0;
            let flipRatio = 0.85;
            let frameNr = 0;
            let _overRelaxation = 1.9;
            let compensateDrift = true;
            let separateParticles = true;
            let paused = false;

            let overRelaxation = 1.9;
            let dt = 1.0 / 60.0;
            let numPressureIters = 10;
            let numParticleIters = 1;

            let res = 23.0;
            let tankHeight = 1.0 * simHeight;
            let tankWidth = 1.0 * simWidth;
            let h = tankHeight / res;
            let density = 1000.0;

            let relWaterHeight = 0.8;
            let relWaterWidth = 0.6;

            let r = 0.5 * h;
            let dx = 2.0 * r;
            let dy = sqrtf(3.0) / 2.0 * dx;

            let numX = floorf((relWaterWidth * tankWidth - 2.0 * h - 2.0 * r) / dx);
            let numY = floorf((relWaterHeight * tankHeight - 2.0 * h - 2.0 * r) / dy);
            let maxParticles = (numX * numY) as i32;

            let mut fluid = FlipFluid::new(density, tankWidth, tankHeight, h, r, maxParticles);
            fluid.numParticles = particles;

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
            }
        }
        pub fn pause(&mut self) {
            self.paused = true;
        }
        pub fn unpause(&mut self) {
            self.paused = false;
        }
        pub fn is_paused(&self) -> bool {
            self.paused
        }
        pub fn get_num_particles(&self) -> i32 {
            self.fluid.numParticles
        }
        pub fn particle_add(&mut self, add: i32, max: i32) {
            let current = self.fluid.numParticles;
            if current + add < 0 {
                self.fluid.numParticles = 0;
            } else if current + add > max {
                self.fluid.numParticles = max;
            } else {
                self.fluid.numParticles = current + add;
            }
        }
        pub fn simulate(&mut self) {
            self.fluid.cellType.fill(CellType::AIR_CELL);
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
        pub fn set_num_particles(&mut self, num_particles: i32) {
            self.fluid.numParticles = num_particles;
        }
        pub fn set_gravity(&mut self, accel_measurment: [f32; 2]) {
            self.xGravity = accel_measurment[0];
            self.yGravity = accel_measurment[1];
        }
        pub fn get_output(&mut self) -> [[bool; 21]; 21] {
            let mut output_frame: [[bool; 21]; 21] = [[false; 21]; 21];
            for i in 1..22 {
                for j in 1..22 {
                    if self.fluid.cellType[i * 23 + j] == CellType::FLUID_CELL {
                        output_frame[i - 1][j - 1] = true;
                    }
                }
            }
            output_frame
        }
    }

    #[unsafe(no_mangle)]
    pub extern "C" fn setup_scene_wrapper(num_particles: i32) -> *mut Scene {
        let scene = Box::new(Scene::setupScene(num_particles));
        Box::into_raw(scene)
    }

    #[unsafe(no_mangle)]
    pub extern "C" fn destroy_scene_wrapper(ptr: *mut Scene) {
        if ptr.is_null() {
            return;
        }
        unsafe {
            drop(Box::from_raw(ptr));
        }
    }

    #[unsafe(no_mangle)]
    pub extern "C" fn simulate_wrapper(ptr: *mut Scene) {
        if ptr.is_null() {
            return;
        }
        let scene = unsafe { &mut *ptr };
        scene.simulate();
    }

    #[unsafe(no_mangle)]
    pub extern "C" fn set_gravity_wrapper(ptr: *mut Scene, x_gravity: f32, y_gravity: f32) {
        if ptr.is_null() {
            return;
        }
        let scene = unsafe { &mut *ptr };
        scene.set_gravity([x_gravity, y_gravity]);
    }

    #[unsafe(no_mangle)]
    pub extern "C" fn get_output_wrapper(ptr: *mut Scene) -> *const bool {
        if ptr.is_null() {
            return ptr::null();
        }
        let scene = unsafe { &mut *ptr };
        let output = scene.get_output();
        let flat_output: Vec<bool> = output.iter().flat_map(|row| row.iter()).cloned().collect();
        let boxed_output = flat_output.into_boxed_slice();
        unsafe { Box::into_raw(boxed_output) as *const bool }
    }

    #[unsafe(no_mangle)]
    pub extern "C" fn get_output_rows() -> i32 {
        21
    }

    #[unsafe(no_mangle)]
    pub extern "C" fn get_output_cols() -> i32 {
        21
    }

    fn clamp<T: PartialOrd>(x: T, min: T, max: T) -> T {
        if x < min {
            min
        } else if x > max {
            max
        } else {
            x
        }
    }
}
