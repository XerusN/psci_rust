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
    // c : f64,
    // cfl : f64,
    // t_f : f64,
    // x_0 : f64,
    // delta : f64,
    // bc : BoundaryCondition,
    frame : usize,
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
    
    let mut file_line = file_lines[6].split(' ');
    let l_x: f64 = file_line.next().unwrap().parse().unwrap();
    let l_y: f64 = file_line.next().unwrap().parse().unwrap();
    
    let frame : usize = file_lines[9].parse().unwrap();
    
    InitInfo {
        n_x : n_x,
        n_y : n_y,
        l_x : l_x,
        l_y : l_y,
        frame : frame,
    }
}

fn write_output(iteration : usize, info : &InitInfo,space_grid : &SpaceGrid, p : &Vec<Vec<f64>>) {
    
    let filename = "target/output_rust/resTECPLOT_".to_owned() + &iteration.to_string() + ".dat";
    let err = fs::create_dir_all("target/output_rust");
    if let Err(err) = err  {
        println!("Error while creating the output folder: {:?}", err);
    }

    let mut file = fs::File::create(filename).unwrap();
    
    let _ = writeln!(file, "TITLE = \"ETAPE 5 Rust\"");
    let _ = writeln!(file, "VARIABLES = \"X\" \"Y\" \"P\"");
    let _ = writeln!(file, "ZONE T= \"0   seconds\", I={0:?}, J={1:?}, DATAPACKING=POINT", info.n_x, info.n_y);
    
    for i in 0..info.n_x {
        
        for j in 0..info.n_y {
            
            let _= writeln!(file, "{0} {1} {2}", fmt_f64(&space_grid.x[i], &25, &16, &3), fmt_f64(&space_grid.y[j], &25, &16, &3), fmt_f64(&p[i][j], &25, &16, &3));
            
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
            
            if i == 0 || i == info.n_x-1 || j == 0 || j == info.n_y-1 {
                
                a[k][k] = 1.0
                
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

fn init_p(info : &InitInfo, space_grid : &SpaceGrid) -> Vec<Vec<f64>> {
    
    let mut p : Vec<Vec<f64>> = Vec::new();
    
    for _i in 0..info.n_x {
        let mut p_line : Vec<f64> = Vec::with_capacity(info.n_y);
        for _j in 0..info.n_y {
            p_line.push(0.0)
        }
        p.push(p_line);
    };
    
    
    
    for i in 0..info.n_y {
        p[0][i] = 1.0;
        p[info.n_x-1][i] = space_grid.y[i].exp();
    }
    for i in 0..info.n_x {
        p[i][0] = 1.0 - 2.0*space_grid.x[i].sinh();
        p[i][info.n_y-1] = space_grid.x[i].exp() - 2.0*space_grid.x[i].sinh();
    }
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



fn init_b(info : &InitInfo, space_grid : &SpaceGrid) -> Vec<f64> {
    
    let k_max = info.n_x*info.n_y;
    let mut b = Vec::with_capacity(k_max);
    
    for j in 0..info.n_y {
        for i in 0..info.n_x {
            if i == 0 {
                b.push(1.0);
            } else if i == info.n_x-1 {
                b.push(space_grid.y[j].exp() - 2.0*(1.0 as f64).sinh());
            }else if j == 0 {
                b.push(1.0 - 2.0*space_grid.x[i].sinh());
            } else if j == info.n_y-1 {
                b.push(space_grid.x[i].exp() - 2.0*space_grid.x[i].sinh())
            } else {
                b.push((space_grid.x[i].powf(2.0) + space_grid.y[j].powf(2.0))*(space_grid.x[i]*space_grid.y[j]).exp() - 2.0*space_grid.x[i].sinh())   
            }   
        }   
    }   
    b
}



fn norm_2(vector : &Vec<f64>) -> f64 {
    let n = vector.len();
    
    let mut norm : f64 = 0.0;
    
    for i in 0..n {       
        norm += vector[i].powf(2.0);
    }
    
    norm.sqrt()
}

fn jacobi_method(info : &InitInfo, a : &Vec<Vec<f64>>, b : &Vec<f64>, p_vec : &mut Vec<f64>, r : &mut Vec<f64>, space_grid : &SpaceGrid, p : &mut Vec<Vec<f64>>) {
    let r_tol = 1e-5;
    
    let k_max = p_vec.len();
    
    p_vec.fill(0.0);
    
    for i in 0..k_max {
        r[i] = -b[i];
        for j in 0..k_max {
            r[i] += a[i][j]*p_vec[j];
        }
    }

    let r_norm_0 = norm_2(r);
    let mut r_norm = r_norm_0;
    let mut iteration = 0;
    
    
    println!("{:?}", r_norm_0);
    
    // println!{"{}", r_norm_0};
    while r_norm > r_norm_0 * r_tol {
        for i in 0..k_max {
            p_vec[i] -= r[i] / a[i][i];
        }
        
        for i in 0..k_max {
            r[i] = -b[i];
            for j in 0..k_max {
                r[i] += a[i][j] * p_vec[j];
            }
        }
        
        r_norm = norm_2(r);
        iteration += 1;
        // if iteration % info.frame == 0 {
        //     matrixize_vector(&info, p, &p_vec);
        //     write_output(iteration, info, space_grid, p);
        // }
        println!("{:?}", iteration);
    }
    
    
    matrixize_vector(&info, p, &p_vec);
    write_output(iteration, info, space_grid, p);    
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
    let b = init_b(&info, &space_grid);
    
    let mut p = init_p(&info, &space_grid);
    write_output(0, &info, &space_grid, &p);
    let mut p_vec : Vec<f64> = Vec::with_capacity(info.n_x*info.n_y);
    let mut r : Vec<f64> = Vec::with_capacity(info.n_x*info.n_y);
    for _i in 0..(info.n_x*info.n_y) {p_vec.push(0.0); r.push(0.0);}
    
    vectorize_matrix(&info, &p, &mut p_vec);
    jacobi_method(&info, &a, &b, &mut p_vec, &mut r, &space_grid, &mut p);
    
}
