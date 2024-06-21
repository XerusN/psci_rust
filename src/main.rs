use core::panic;
use std::fs;
use std::io::{BufRead, BufReader};
use std::io::prelude::*;

struct SpaceGrid {
    x : Vec<f64>,
    y : Vec<f64>,
}

struct InitInfo {
    n_x : usize,
    n_y : usize,
    l_x : f64,
    l_y : f64,
    viscosity : f64,
    density : f64,
    cfl : f64,
    fo : f64,
    re : f64,
    t_f : f64,
    frame : usize,
    scheme : usize,
}

fn fmt_f64(num: &f64, width: &usize, precision: &usize, exp_pad: &usize) -> String {
    let mut num = format!("{:.precision$e}", num, precision = precision);
    // Safe to `unwrap` as `num` is guaranteed to contain `'e'`
    let exp = num.split_off(num.find('e').unwrap());

    let (sign, exp) = if exp.starts_with("e-") {
        ('-', &exp[2..])
    } else {
        ('+', &exp[1..])
    };
    num.push_str(&format!("e{}{:0>pad$}", sign, exp, pad = exp_pad));

    format!("{:>width$}", num, width = width)
}

fn read_info(file_path : &str) -> InitInfo{
    
    let f = fs::File::open(file_path).unwrap();
    let reader = BufReader::new(f);
    let file_lines = reader.lines().into_iter().flatten().collect::<Vec<_>>();
    
    let mut file_line = file_lines[3].split(' ');
    let n_x: usize = file_line.next().unwrap().parse().unwrap();
    let n_y: usize = file_line.next().unwrap().parse().unwrap();
    
    let t_f : f64 = file_lines[6].parse().unwrap();
    
    let mut file_line = file_lines[9].split(' ');
    let l_x: f64 = file_line.next().unwrap().parse().unwrap();
    let l_y: f64 = file_line.next().unwrap().parse().unwrap();
    
    let frame : usize = file_lines[12].parse().unwrap();
    
    let density : f64 = file_lines[15].parse().unwrap();
    
    let cfl : f64 = file_lines[17].parse().unwrap();
    
    let fo : f64 = file_lines[19].parse().unwrap();
    
    let re : f64 = file_lines[22].parse().unwrap();
    
    let scheme : usize = file_lines[25].parse().unwrap();
    
    let viscosity : f64 = 1.0*1.0/re;
    
    InitInfo {
        n_x : n_x,
        n_y : n_y,
        l_x : l_x,
        l_y : l_y,
        viscosity : viscosity,
        density : density,
        cfl : cfl,
        fo : fo,
        re : re,
        t_f : t_f,
        frame : frame,
        scheme : scheme,
    }
    
}

fn write_output(iteration : usize, info : &InitInfo,space_grid : &SpaceGrid, u : &Vec<Vec<f64>>, v : &Vec<Vec<f64>>, p : &Vec<Vec<f64>>) {
    
    let filename = "target/output_rust/resTECPLOT_".to_owned() + &iteration.to_string() + ".dat";
    let err = fs::create_dir_all("target/output_rust");
    if let Err(err) = err  {
        println!("Error while creating the output folder: {:?}", err);
    }

    let mut file = fs::File::create(filename).unwrap();
    
    let _ = writeln!(file, "TITLE = \"ETAPE 5 Rust\"");
    let _ = writeln!(file, "VARIABLES = \"X\" \"Y\" \"U\" \"V\" \"P\"");
    let _ = writeln!(file, "ZONE T= \"0   seconds\", I={0:?}, J={1:?}, DATAPACKING=POINT", info.n_x, info.n_y);
    
    for i in 0..info.n_x {
        
        for j in 0..info.n_y {
            
            let _= writeln!(file, "{0} {1} {2} {3} {4}", fmt_f64(&space_grid.x[i], &25, &16, &3), fmt_f64(&space_grid.y[j], &25, &16, &3), fmt_f64(&u[i][j], &25, &16, &3), fmt_f64(&v[i][j], &25, &16, &3), fmt_f64(&p[i][j], &25, &16, &3));
            
        }
        
    }
    
}



fn gen_space_grid(info : &InitInfo) -> SpaceGrid {
    
    let dx : f64 = info.l_x/((info.n_x as f64) - 1.0);
    let dy : f64 = info.l_y/((info.n_y as f64) - 1.0);
    
    let mut space_grid : SpaceGrid = SpaceGrid {
        x : Vec::new(),
        y : Vec::new(),
    };
    
    
    for i in 0..info.n_x {
        space_grid.x.push(dx*(i as f64));
    }
    for i in 0..info.n_y {
        space_grid.y.push(dy*(i as f64));
    }
    space_grid
}



fn init_a(info : &InitInfo) -> Vec<Vec<f64>> {
    
    let k_max = info.n_x*info.n_y;
    let mut a = Vec::new();
    
    for _i in 0..k_max {
        let mut a_line : Vec<f64> = Vec::with_capacity(k_max);
        for _j in 0..k_max {
            a_line.push(0.0);
        }
        a.push(a_line);
    };
    
    let dx_2_inv = 1.0/(info.l_x/((info.n_x as f64) - 1.0)).powf(2.0);
    let dy_2_inv = 1.0/(info.l_y/((info.n_y as f64) - 1.0)).powf(2.0);
    
    for j in 0..info.n_y {
        
        for i in 0..info.n_x {
            
            let k = j*info.n_x + i;
            
            if i == 0 {
                a[k][k] = 1.0;
                a[k][k+1] = -1.0;
            } else if i == info.n_x-1 {
                a[k][k] = 1.0;
                a[k][k-1] = -1.0;
            }else if j == 0 {
                a[k][k] = 1.0;
                a[k][k+info.n_x] = -1.0;
            } else if j == info.n_y-1 {
                a[k][k] = 1.0;
                a[k][k-info.n_x] = -1.0;
            } else {
                a[k][k] = -2.0*(dx_2_inv + dy_2_inv);
                a[k][k-1] = dx_2_inv;
                a[k][k+1] = dx_2_inv;
                a[k][k-info.n_x] = dy_2_inv;
                a[k][k+info.n_x] = dy_2_inv;
            }
            
        }
        
    }
    
    a
}



fn init_u(info : &InitInfo) -> Vec<Vec<f64>> {
    
    let mut u : Vec<Vec<f64>> = Vec::new();
    
    for _i in 0..info.n_x {
        let mut u_line : Vec<f64> = Vec::with_capacity(info.n_y);
        for _j in 0..(info.n_y - 1) {
            u_line.push(0.0)
        }
        u_line.push(1.0);
        u.push(u_line);
    };
    
    u
}



fn init_v(info : &InitInfo) -> Vec<Vec<f64>> {
    
    let mut v : Vec<Vec<f64>> = Vec::new();
    
    for _i in 0..info.n_x {
        let mut v_line : Vec<f64> = Vec::with_capacity(info.n_y);
        for _j in 0..info.n_y {
            v_line.push(0.0)
        }
        v.push(v_line);
    };
    
    v
}

fn init_p(info : &InitInfo) -> Vec<Vec<f64>> {
    
    let mut p : Vec<Vec<f64>> = Vec::new();
    
    for _i in 0..info.n_x {
        let mut p_line : Vec<f64> = Vec::with_capacity(info.n_y);
        for _j in 0..info.n_y {
            p_line.push(0.0)
        }
        p.push(p_line);
    };
    
    p
}



fn vectorize_matrix(info : &InitInfo, matrix : &Vec<Vec<f64>>, vector : &mut Vec<f64>) {

    for i in 0..info.n_x {
        
        for j in 0..info.n_y {
            
            let k = j*info.n_x + i;
            
            vector[k] = matrix[i][j]
            
        }   
    }   
}



fn matrixize_vector(info : &InitInfo, matrix : &mut Vec<Vec<f64>>, vector : &Vec<f64>) {
    
    for i in 0..info.n_x {
        
        for j in 0..info.n_y {
            
            let k = j*info.n_x + i;
            
            matrix[i][j] = vector[k]
            
        }
        
    }
    
}


fn norm_2(vector : &Vec<f64>) -> f64 {
    
    let mut norm : f64 = 0.0;
    
    for i in 0..vector.len() {       
        norm += vector[i].powf(2.0);
    }
    
    norm.sqrt()
}

fn compute_pressure_integral(info : &InitInfo, p_vec : &mut Vec<f64>) -> f64 {
    
    let mut integral : f64 = 0.0;
    
    let dx : f64 = info.l_x/((info.n_x as f64) - 1.0);
    let dy : f64 = info.l_y/((info.n_y as f64) - 1.0);
    
    for i in 0..info.n_x {
        for j in 0..info.n_y {
            
            let k = j*info.n_x + i;
            
            if ((i == 0) | (i == info.n_x-1)) ^ ((j == 0) | (j == info.n_y-1)) {
                integral += p_vec[k]*dx*dy*0.25;
            }
            else if (i == 0) | (i == info.n_x-1) | (j == 0) | (j == info.n_y-1) {
                integral += p_vec[k]*dx*dy*0.5;
            }
            else {
                integral += p_vec[k]*dx*dy;
            }
        }
    }
    
    integral /= info.l_x*info.l_y;
    
    for i in 0..info.n_x {
        for j in 0..info.n_y {
            
            let k = j*info.n_x + i;
            
            p_vec[k] -= integral;
        }
    }
    
    integral
    
}

fn jacobi_method(info : &InitInfo, a : &Vec<Vec<f64>>, b : &Vec<f64>, p_vec : &mut Vec<f64>, p_vec_temp : &mut Vec<f64>, r : &mut Vec<f64>) {
    
    let epsilon = 1e-3;
    let n_max: usize = 20000;
    
    let k_max = p_vec.len();
    
    let mut integral : f64 = compute_pressure_integral(info, p_vec);
    
    for i in 0..k_max {
        r[i] = -b[i];
        for j in 0..k_max {
            r[i] = r[i] + a[i][j]*p_vec[j];
        }
    }
    
    let mut upper_norm : f64 = norm_2(&r);
    let mut lower_norm : f64 = norm_2(&b);
    
    let mut iteration = 0;
    
    
    while upper_norm > epsilon*lower_norm {
        
        for i in 0..k_max {
            p_vec_temp[i] = p_vec[i];
        }
        
        for i in 0..k_max {
            p_vec[i] = p_vec_temp[i] - r[i] / a[i][i];
        }
        
        integral = compute_pressure_integral(info, p_vec);
        
        for i in 0..k_max {
            r[i] = -b[i];
            for j in 0..k_max {
                r[i] = r[i] + a[i][j]*p_vec[j];
            }
            p_vec_temp[i] = p_vec[i] - p_vec_temp[i]
        }
        
        iteration += 1;
        
        upper_norm = norm_2(&p_vec_temp);
        lower_norm = norm_2(&p_vec);
        //upper_norm = norm_2(&r);
        //lower_norm = norm_2(&b);
        
        if iteration >= n_max {
            panic!("Nombre d'iteration max atteint sur jacobi : {:?}", n_max);
        }
        
        if iteration % 100 == 0 {
            println!("convergence : {:?} | integral : {:?}", upper_norm/lower_norm, integral);
        }
    }
    println!("jacobi : {:?} iterations | integral = {:?}", iteration, integral);
    
}



fn compute_b(info : &InitInfo, dt : f64, b : &mut Vec<f64>, u_temp : &Vec<Vec<f64>>, v_temp : &Vec<Vec<f64>>){
    
    let mut k: usize;
    
    let dx : f64 = info.l_x/((info.n_x as f64) - 1.0);
    let dy : f64 = info.l_y/((info.n_y as f64) - 1.0);
    
    for j in 0..info.n_y {
        for i in 0..info.n_x {
            k = info.n_x * j + i;
            if i == 0 {
                b[k] = 0.0;
            } else if i == info.n_x-1 {
                b[k] = 0.0;
            }else if j == 0 {
                b[k] = 0.0;
            } else if j == info.n_y-1 {
                b[k] = 0.0;
            } else {
                b[k] = (u_temp[i+1][j] - u_temp[i-1][j])/(2.0*dx) + (v_temp[i][j+1] - v_temp[i][j-1])/(2.0*dy);
                b[k] *= info.density/dt;
            }   
        }
    }
}



fn compute_time_step(info : &InitInfo, u : &Vec<Vec<f64>>, v : &Vec<Vec<f64>>) -> f64 {
    
    let dx : f64 = info.l_x/((info.n_x as f64) - 1.0);
    let dy : f64 = info.l_y/((info.n_y as f64) - 1.0);
    
    let mut dt : f64 = info.fo*dx.powf(2.0)/info.viscosity;
    
    for i in 0..info.n_x {
        for j in 0..info.n_y {
            if u[i][j].abs() > 0.0 {
                if (info.cfl*dx/u[i][j]).abs() < dt {
                    dt = (info.cfl*dx/u[i][j]).abs()
                }
            }
            if v[i][j].abs() > 0.0 {
                if (info.cfl*dy/v[i][j]).abs() < dt {
                    dt = (info.cfl*dy/v[i][j]).abs()
                }
            }
        }
        
    }
    
    dt
}



fn speed_guess(info : &InitInfo, u : &Vec<Vec<f64>>, u_temp : &mut Vec<Vec<f64>>, v : &Vec<Vec<f64>>, v_temp : &mut Vec<Vec<f64>>, dt : f64) {
    
    let dx : f64 = info.l_x/((info.n_x as f64) - 1.0);
    let dy : f64 = info.l_y/((info.n_y as f64) - 1.0);
    
    u_temp[0].fill(0.0);
    u_temp[info.n_x-1].fill(0.0);
    v_temp[0].fill(0.0);
    v_temp[info.n_x-1].fill(0.0);
    for i in 0..info.n_x {
        u_temp[i][0] = 0.0;
        u_temp[i][info.n_y-1] = 1.0;
        v_temp[i][0] = 0.0;
        v_temp[i][info.n_y-1] = 0.0;
    }
    
    
    match info.scheme {
        
        0 => {
            let mut upwind_x : f64;
            let mut upwind_y : f64;
            
            for i in 1..info.n_x-1 {
                for j in 1..info.n_y-1 {
                    
                    if u[i][j] > 0.0 {
                        upwind_x = 1.0
                    } else {
                        upwind_x = 0.0
                    }
                    
                    if v[i][j] > 0.0 {
                        upwind_y = 1.0
                    } else {
                        upwind_y = 0.0
                    }
                    
                    u_temp[i][j] = u[i][j]*(1.0 - 2.0*info.viscosity*dt*(1.0/dx.powf(2.0) + 1.0/dy.powf(2.0)) - u[i][j]*dt/dx*2.0*(upwind_x - 0.5) - v[i][j]*dt/dy*2.0*(upwind_y - 0.5)) + u[i-1][j]*dt*(info.viscosity/dx.powf(2.0) + u[i][j]/dx*upwind_x) + u[i][j-1]*dt*(info.viscosity/dy.powf(2.0) + v[i][j]/dy*upwind_y) + u[i+1][j]*dt*(info.viscosity/dx.powf(2.0) + u[i][j]/dx*(1.0 - upwind_x)) + u[i][j+1]*dt*(info.viscosity/dy.powf(2.0) + v[i][j]/dy*(1.0 - upwind_y));
                    v_temp[i][j] = v[i][j]*(1.0 - 2.0*info.viscosity*dt*(1.0/dx.powf(2.0) + 1.0/dy.powf(2.0)) - u[i][j]*dt/dx*2.0*(upwind_x - 0.5) - v[i][j]*dt/dy*2.0*(upwind_y - 0.5)) + v[i-1][j]*dt*(info.viscosity/dx.powf(2.0) + u[i][j]/dx*upwind_x) + v[i][j-1]*dt*(info.viscosity/dy.powf(2.0) + v[i][j]/dy*upwind_y) + v[i+1][j]*dt*(info.viscosity/dx.powf(2.0) + u[i][j]/dx*(1.0 - upwind_x)) + v[i][j+1]*dt*(info.viscosity/dy.powf(2.0) + v[i][j]/dy*(1.0 - upwind_y));
                    
                    // if (u_temp[i][j] > 1.0) | (v_temp[i][j] > 1.0) {
                    //     panic!("speed guess")
                    // }
                    
                }
            
            }
        },
        
        _ => {
            
            for i in 1..info.n_x-1 {
                for j in 1..info.n_y-1 {
                    
                    u_temp[i][j] = u[i][j]*(1.0 - 2.0*info.viscosity*dt*(1.0/dx.powf(2.0) + 1.0/dy.powf(2.0))) + u[i-1][j]*dt*(info.viscosity/dx.powf(2.0) + u[i][j]/(2.0*dx)) + u[i][j-1]*dt*(info.viscosity/dy.powf(2.0) + v[i][j]/(2.0*dy)) + u[i+1][j]*dt*(info.viscosity/dx.powf(2.0) - u[i][j]/(2.0*dx)) + u[i][j+1]*dt*(info.viscosity/dy.powf(2.0) - v[i][j]/(2.0*dy));
                    v_temp[i][j] = v[i][j]*(1.0 - 2.0*info.viscosity*dt*(1.0/dx.powf(2.0) + 1.0/dy.powf(2.0))) + v[i-1][j]*dt*(info.viscosity/dx.powf(2.0) + u[i][j]/(2.0*dx)) + v[i][j-1]*dt*(info.viscosity/dy.powf(2.0) + v[i][j]/(2.0*dy)) + v[i+1][j]*dt*(info.viscosity/dx.powf(2.0) - u[i][j]/(2.0*dx)) + v[i][j+1]*dt*(info.viscosity/dy.powf(2.0) - v[i][j]/(2.0*dy));

                }
            
            }
        },
        
        
    };
    
    
}



fn compute_pressure(info : &InitInfo, dt : f64, p : &mut Vec<Vec<f64>>, u_temp : &Vec<Vec<f64>>, v_temp : &Vec<Vec<f64>>, a : &Vec<Vec<f64>>, b : &mut Vec<f64>, p_vec : &mut Vec<f64>, p_vec_temp : &mut Vec<f64>, r : &mut Vec<f64>) {
    
    compute_b(info, dt, b, u_temp, v_temp);
    
    vectorize_matrix(info, p, p_vec);
    
    jacobi_method(info, a, b, p_vec, p_vec_temp, r);
    
    matrixize_vector(info, p, p_vec);
    
}



fn adjust_speed(info : &InitInfo, dt : f64, p : &Vec<Vec<f64>>, u : &mut Vec<Vec<f64>>, u_temp : &Vec<Vec<f64>>, v : &mut Vec<Vec<f64>>, v_temp : &Vec<Vec<f64>>) {
    
    let dx : f64 = info.l_x/((info.n_x as f64) - 1.0);
    let dy : f64 = info.l_y/((info.n_y as f64) - 1.0);
    
    u[0].fill(0.0);
    u[info.n_x-1].fill(0.0);
    v[0].fill(0.0);
    v[info.n_x-1].fill(0.0);
    for i in 0..info.n_x {
        u[i][0] = 0.0;
        u[i][info.n_y-1] = 1.0;
        v[i][0] = 0.0;
        v[i][info.n_y-1] = 0.0;
    }
    
    for i in 1..info.n_x-1 {
        for j in 1..info.n_y-1 {
            
            u[i][j] = u_temp[i][j] - dt/info.density*(p[i+1][j] - p[i-1][j])/(2.0*dx);
            v[i][j] = v_temp[i][j] - dt/info.density*(p[i][j+1] - p[i][j-1])/(2.0*dy);
            
            // if (u[i][j] > 1.0) | (v[i][j] > 1.0) {
            //     panic!("adjust speed")
            // }
        }
    }
}

fn main() {
    
    let e = fs::remove_dir_all("./target/output_rust/");
    if let Err(e) = e {
        println!("Error while removing files from output_rust : {:?}", e);
    }

    let info : InitInfo = read_info("init/input.dat");
    
    let space_grid = gen_space_grid(&info);
    // for i in space_grid {
    //     println!("{}", i.x);
    // }
    let a = init_a(&info);
    let mut b = Vec::with_capacity(info.n_x*info.n_y);
    
    
    
    let mut u = init_u(&info);
    let mut u_temp = init_u(&info);
    let mut v = init_v(&info);
    let mut v_temp = init_v(&info);
    let mut p = init_p(&info);
    let mut p_vec : Vec<f64> = Vec::with_capacity(info.n_x*info.n_y);
    let mut p_vec_temp : Vec<f64> = Vec::with_capacity(info.n_x*info.n_y);
    let mut r : Vec<f64> = Vec::with_capacity(info.n_x*info.n_y);
    for _i in 0..(info.n_x*info.n_y) {p_vec.push(0.0); r.push(0.0); b.push(0.0); p_vec_temp.push(0.0);}
    
    let mut t: f64 = 0.0;
    
    let mut dt: f64;
    
    let mut last_iteration = false;
    let mut iteration = 0;
    
    write_output(iteration, &info, &space_grid, &u, &v, &p);
    
    while last_iteration == false {
        
        dt = compute_time_step(&info, &u, &v);
        
        if t + dt >= info.t_f {
            last_iteration = true;
            break;
        }
        println!("dt = {:?}", dt);
        
        speed_guess(&info, &u, &mut u_temp, &v, &mut v_temp, dt);
        
        //write_output(3, &info, &space_grid, &u_temp, &v_temp, &p);
        // if iteration % info.frame == 0 {
        //     write_output(1, &info, &space_grid, &u_temp, &v_temp  , &p);
        // }
        
        // println!("OK");
        compute_pressure(&info, dt, &mut p, &u_temp, &v_temp, &a, &mut b, &mut p_vec, &mut p_vec_temp, &mut r);
        
        adjust_speed(&info, dt, &p, &mut u, &u_temp, &mut v, &v_temp);
        
        t += dt;
        iteration += 1;
        
        println!("t = {:?}", t);
        
        println!("-------------");
        
        if iteration % info.frame == 0 {
            write_output(iteration, &info, &space_grid, &u, &v, &p);
        }
    }
    
    write_output(iteration, &info, &space_grid, &u, &v, &p);
    
}
