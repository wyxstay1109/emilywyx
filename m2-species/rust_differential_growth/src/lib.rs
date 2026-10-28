// =========================================================================
// 白蜃蠕蛇 (Gypsochiton serpentinus) - Rust WebAssembly 高性能運算模組
// 負責：
// 1. 微分生長 (Differential Growth) 高階生理微血管網絡物理模擬
// 2. 血管新生 (Angiogenesis) 與微血管床網狀吻合 (Capillary Anastomosis) 運算
// 3. 38節脊椎反向運動學 (Inverse Kinematics) 幾何流暢傳播
// =========================================================================

use std::cell::RefCell;

thread_local! {
    static RESULT_BUFFER: RefCell<Vec<f32>> = RefCell::new(Vec::new());
    static IK_OUTPUT: RefCell<Vec<f32>> = RefCell::new(Vec::new());
}

// 輕量級快速 PRNG (Xorshift32)，不依賴外部 Crate
struct SimpleRng {
    state: u32,
}

impl SimpleRng {
    fn new(seed: u32) -> Self {
        Self { state: if seed == 0 { 0x12345678 } else { seed } }
    }
    fn next_u32(&mut self) -> u32 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 17;
        x ^= x << 5;
        self.state = x;
        x
    }
    fn next_f32(&mut self) -> f32 {
        (self.next_u32() as f32) / (u32::MAX as f32)
    }
    fn range_f32(&mut self, min: f32, max: f32) -> f32 {
        min + self.next_f32() * (max - min)
    }
}

#[derive(Clone, Copy)]
struct Node {
    x: f32,
    y: f32,
    vx: f32,
    vy: f32,
    order: u8, // 0: 主幹血管, 1: 次級分枝, 2: 細小微血管 (Capillary)
}

struct Path {
    nodes: Vec<Node>,
    order: u8,
}

#[derive(Clone, Copy)]
struct Anastomosis {
    x1: f32,
    y1: f32,
    cpx: f32,
    cpy: f32,
    x2: f32,
    y2: f32,
    is_capillary: u8,
}

// 執行體表微分生長微血管運算 (Differential Growth)
#[no_mangle]
pub extern "C" fn generate_differential_growth_vasculature(seed: u32) -> usize {
    let mut rng = SimpleRng::new(seed);

    let width: f32 = 2048.0;
    let height: f32 = 1024.0;

    let max_edge: f32 = 24.0;
    let target_dist: f32 = 14.0;
    let repel_radius: f32 = 30.0;
    let repel_force: f32 = 0.95;
    let spring_force: f32 = 0.48;
    let smooth_force: f32 = 0.38;

    let mut paths: Vec<Path> = Vec::with_capacity(72);

    // 1. 生成遍布體表的主幹與側腹種子路徑 (共 22 條縱向分佈主脈，徹底覆蓋體表四周)
    let num_seeds: usize = 22;
    for p in 0..num_seeds {
        let mut nodes = Vec::with_capacity(36);
        let base_y = (p as f32 + 0.5) * (height / num_seeds as f32);
        let num_nodes = 24;
        let order = if p % 2 == 0 { 0 } else { 1 };

        for n in 0..num_nodes {
            let x = (n as f32 / (num_nodes - 1) as f32) * (width - 60.0) + 30.0;
            let y = base_y + rng.range_f32(-16.0, 16.0);
            nodes.push(Node {
                x,
                y,
                vx: 0.0,
                vy: 0.0,
                order,
            });
        }
        paths.push(Path {
            nodes,
            order,
        });
    }

    // 2. 空間加速哈希網格 (Spatial Hashing Grid)
    let cell_size: f32 = 34.0;
    let grid_cols = (width / cell_size).ceil() as usize + 2;
    let grid_rows = (height / cell_size).ceil() as usize + 2;
    let mut grid: Vec<Vec<(usize, usize)>> = vec![Vec::new(); grid_cols * grid_rows];

    // 3. 微分生長迭代模擬 (Differential Growth Physics Simulation - 80 步精密演算)
    let total_steps = 80;
    for step in 0..total_steps {
        // (A) 節點細分 (Subdivision) - 沿血管生長延展
        let current_path_count = paths.len();
        for pi in 0..current_path_count {
            let mut i = 0;
            while i < paths[pi].nodes.len().saturating_sub(1) {
                let dx = paths[pi].nodes[i + 1].x - paths[pi].nodes[i].x;
                let dy = paths[pi].nodes[i + 1].y - paths[pi].nodes[i].y;
                let dist = (dx * dx + dy * dy).sqrt();

                if dist > max_edge {
                    let mid_x = (paths[pi].nodes[i].x + paths[pi].nodes[i + 1].x) * 0.5 + rng.range_f32(-2.2, 2.2);
                    let mid_y = (paths[pi].nodes[i].y + paths[pi].nodes[i + 1].y) * 0.5 + rng.range_f32(-2.2, 2.2);
                    let order = paths[pi].nodes[i].order;

                    paths[pi].nodes.insert(i + 1, Node {
                        x: mid_x,
                        y: mid_y,
                        vx: 0.0,
                        vy: 0.0,
                        order,
                    });
                    i += 1;
                }
                i += 1;
            }

            // (B) 血管新生 (Angiogenic Sprouting): 滋生鮮紅細小微血管 (Capillaries)
            if step > 5 && step % 5 == 0 && paths.len() < 58 {
                let n_len = paths[pi].nodes.len();
                if n_len > 6 {
                    let branch_idx = (rng.next_u32() as usize % (n_len - 4)) + 2;
                    let b_node = paths[pi].nodes[branch_idx];
                    let sign = if rng.next_f32() > 0.5 { 1.0 } else { -1.0 };
                    let angle = sign * (std::f32::consts::PI * 0.45) + rng.range_f32(-0.35, 0.35);

                    let b1 = Node {
                        x: b_node.x,
                        y: b_node.y,
                        vx: 0.0,
                        vy: 0.0,
                        order: 2, // 細小微血管
                    };
                    let b2 = Node {
                        x: b_node.x + angle.cos() * target_dist * 0.95,
                        y: b_node.y + angle.sin() * target_dist * 0.95,
                        vx: 0.0,
                        vy: 0.0,
                        order: 2,
                    };
                    paths.push(Path {
                        nodes: vec![b1, b2],
                        order: 2,
                    });
                }
            }
        }

        // (C) 重建空間網格
        for cell in grid.iter_mut() {
            cell.clear();
        }
        for (pi, path) in paths.iter().enumerate() {
            for (ni, node) in path.nodes.iter().enumerate() {
                let gx = (node.x / cell_size).max(0.0).min((grid_cols - 1) as f32) as usize;
                let gy = (node.y / cell_size).max(0.0).min((grid_rows - 1) as f32) as usize;
                let cell_idx = gy * grid_cols + gx;
                if cell_idx < grid.len() {
                    grid[cell_idx].push((pi, ni));
                }
            }
        }

        // (D) 微分生長多體力學演算：排斥力、彈簧張力與曲率平滑
        for pi in 0..paths.len() {
            let path_len = paths[pi].nodes.len();
            for ni in 0..path_len {
                let node = paths[pi].nodes[ni];
                let mut fx: f32 = 0.0;
                let mut fy: f32 = 0.0;

                let gx = (node.x / cell_size).max(0.0).min((grid_cols - 1) as f32) as isize;
                let gy = (node.y / cell_size).max(0.0).min((grid_rows - 1) as f32) as isize;

                // 排斥力 (Non-overlapping Self-avoiding Walk)
                for ox in -1..=1 {
                    for oy in -1..=1 {
                        let nx = gx + ox;
                        let ny = gy + oy;
                        if nx >= 0 && nx < grid_cols as isize && ny >= 0 && ny < grid_rows as isize {
                            let c_idx = (ny as usize) * grid_cols + (nx as usize);
                            for &(other_pi, other_ni) in &grid[c_idx] {
                                if other_pi == pi && (other_ni as isize - ni as isize).abs() <= 1 {
                                    continue;
                                }
                                let other = paths[other_pi].nodes[other_ni];
                                let dx = node.x - other.x;
                                let dy = node.y - other.y;
                                let d = (dx * dx + dy * dy).sqrt();
                                if d > 0.001 && d < repel_radius {
                                    let force = (1.0 - d / repel_radius) * repel_force;
                                    fx += (dx / d) * force;
                                    fy += (dy / d) * force;
                                }
                            }
                        }
                    }
                }

                // 彈簧連續性張力 (Spring tension)
                if ni > 0 {
                    let prev = paths[pi].nodes[ni - 1];
                    let dx = node.x - prev.x;
                    let dy = node.y - prev.y;
                    let d = (dx * dx + dy * dy).sqrt();
                    if d > 0.001 {
                        let s_force = (d - target_dist) * spring_force;
                        fx -= (dx / d) * s_force;
                        fy -= (dy / d) * s_force;
                    }
                }
                if ni + 1 < path_len {
                    let next = paths[pi].nodes[ni + 1];
                    let dx = node.x - next.x;
                    let dy = node.y - next.y;
                    let d = (dx * dx + dy * dy).sqrt();
                    if d > 0.001 {
                        let s_force = (d - target_dist) * spring_force;
                        fx -= (dx / d) * s_force;
                        fy -= (dy / d) * s_force;
                    }
                }

                // 拉普拉斯曲率平滑 (Laplacian Smoothing)
                if ni > 0 && ni + 1 < path_len {
                    let prev = paths[pi].nodes[ni - 1];
                    let next = paths[pi].nodes[ni + 1];
                    let target_x = (prev.x + next.x) * 0.5;
                    let target_y = (prev.y + next.y) * 0.5;
                    fx += (target_x - node.x) * smooth_force;
                    fy += (target_y - node.y) * smooth_force;
                }

                // 速度阻尼與積分
                let node_ref = &mut paths[pi].nodes[ni];
                node_ref.vx = node_ref.vx * 0.45 + fx * 0.55;
                node_ref.vy = node_ref.vy * 0.45 + fy * 0.55;
            }
        }

        // (E) 位置更新與邊界約束
        for path in paths.iter_mut() {
            for node in path.nodes.iter_mut() {
                node.x = (node.x + node.vx).max(16.0).min(width - 16.0);
                node.y = (node.y + node.vy).max(16.0).min(height - 16.0);
            }
        }
    }

    // 4. 微血管網狀吻合支 (Capillary Anastomosis) - 形成遍布體表肉眼清晰可見的鮮紅微血管網
    let mut anastomoses: Vec<Anastomosis> = Vec::with_capacity(380);
    let num_anast_attempts = 340;
    for _ in 0..num_anast_attempts {
        let p1_idx = rng.next_u32() as usize % paths.len();
        let p2_idx = rng.next_u32() as usize % paths.len();
        if p1_idx == p2_idx || paths[p1_idx].nodes.is_empty() || paths[p2_idx].nodes.is_empty() {
            continue;
        }
        let n1_idx = rng.next_u32() as usize % paths[p1_idx].nodes.len();
        let n2_idx = rng.next_u32() as usize % paths[p2_idx].nodes.len();

        let n1 = paths[p1_idx].nodes[n1_idx];
        let n2 = paths[p2_idx].nodes[n2_idx];

        let dx = n1.x - n2.x;
        let dy = n1.y - n2.y;
        let dist = (dx * dx + dy * dy).sqrt();

        if dist > 12.0 && dist < 90.0 {
            let cpx = (n1.x + n2.x) * 0.5 + rng.range_f32(-14.0, 14.0);
            let cpy = (n1.y + n2.y) * 0.5 + rng.range_f32(-14.0, 14.0);
            let is_capillary = if n1.order == 2 || n2.order == 2 { 1 } else { 0 };

            anastomoses.push(Anastomosis {
                x1: n1.x,
                y1: n1.y,
                cpx,
                cpy,
                x2: n2.x,
                y2: n2.y,
                is_capillary,
            });
        }
    }

    // 5. 序列化至扁平化結果緩衝區 (Flat Binary Layout)
    RESULT_BUFFER.with(|buf_cell| {
        let mut buf = buf_cell.borrow_mut();
        buf.clear();
        // 寫入路徑數量
        buf.push(paths.len() as f32);

        for path in &paths {
            buf.push(path.nodes.len() as f32);
            buf.push(path.order as f32); // 0: 主幹, 1: 分枝, 2: 細小微血管
            for node in &path.nodes {
                buf.push(node.x);
                buf.push(node.y);
            }
        }

        // 寫入吻合微血管數量
        buf.push(anastomoses.len() as f32);
        for an in &anastomoses {
            buf.push(an.x1);
            buf.push(an.y1);
            buf.push(an.cpx);
            buf.push(an.cpy);
            buf.push(an.x2);
            buf.push(an.y2);
            buf.push(an.is_capillary as f32);
        }

        buf.len()
    })
}

// 取得結果緩衝區指標
#[no_mangle]
pub extern "C" fn get_result_buffer_ptr() -> *const f32 {
    RESULT_BUFFER.with(|buf_cell| buf_cell.borrow().as_ptr())
}

// 取得結果緩衝區長度
#[no_mangle]
pub extern "C" fn get_result_buffer_len() -> usize {
    RESULT_BUFFER.with(|buf_cell| buf_cell.borrow().len())
}

// =========================================================================
// 蠕蛇多節反向運動學 (Inverse Kinematics) 與波浪姿態平滑傳播計算
// =========================================================================
#[no_mangle]
pub extern "C" fn compute_worm_ik(
    num_segments: usize,
    segment_len: f32,
    head_x: f32,
    head_y: f32,
    head_z: f32,
    phase: f32,
    amp_deg: f32,
    wave_effect: f32,
    input_coords_ptr: *const f32,
) -> usize {
    let input_slice = unsafe { std::slice::from_raw_parts(input_coords_ptr, num_segments * 3) };

    IK_OUTPUT.with(|out_cell| {
        let mut out = out_cell.borrow_mut();
        out.clear();
        out.reserve(num_segments * 4); // x, y, z, angle

        let amp_rad = amp_deg * (std::f32::consts::PI / 180.0);

        let mut prev_x = head_x;
        let mut prev_y = head_y;
        let mut prev_z = head_z;

        // 第 0 節
        out.push(prev_x);
        out.push(prev_y);
        out.push(prev_z);
        out.push(0.0);

        for i in 1..num_segments {
            let in_x = input_slice[i * 3];
            let in_y = input_slice[i * 3 + 1];
            let in_z = input_slice[i * 3 + 2];

            let dx = in_x - prev_x;
            let dz = in_z - prev_z;
            let natural_angle = dx.atan2(dz);

            let body_wave = (phase - (i as f32) * 0.24).sin() * (amp_rad * 0.55);
            let target_angle = natural_angle + body_wave * wave_effect;

            let cur_x = prev_x + target_angle.sin() * segment_len;
            let cur_z = prev_z + target_angle.cos() * segment_len;
            let cur_y = in_y + (prev_y - in_y) * 0.16;

            let angle = (prev_x - cur_x).atan2(prev_z - cur_z);

            out.push(cur_x);
            out.push(cur_y);
            out.push(cur_z);
            out.push(angle);

            prev_x = cur_x;
            prev_y = cur_y;
            prev_z = cur_z;
        }

        out.len()
    })
}

#[no_mangle]
pub extern "C" fn get_ik_output_ptr() -> *const f32 {
    IK_OUTPUT.with(|out_cell| out_cell.borrow().as_ptr())
}

#[no_mangle]
pub extern "C" fn get_ik_output_len() -> usize {
    IK_OUTPUT.with(|out_cell| out_cell.borrow().len())
}
