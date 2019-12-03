use std::f64;

pub fn trapezoid(start: f64, end: f64, steps: i64, fp: fn(f64) -> f64) -> f64 {
    let h = (end - start) / (steps as f64);

    let mut total = 0.5 * (fp(end) - fp(start));

    for i in 0..steps {
        total += fp(start + h * (i as f64));
    }

    return total * h;
}

pub fn trapezoid_iterative(start: f64, end: f64, _steps: i64, fp: fn(f64) -> f64) -> f64 {
    // adaptive method -- ignore provided 'steps'
    let mut steps = 1;

    let mut old_estimate = 0.0;
    let mut new_estimate = trapezoid_iterative_helper(start, end, steps, fp, old_estimate);

    while (old_estimate - new_estimate) * (old_estimate - new_estimate) > 1e-18 {
        steps += 1;
        old_estimate = new_estimate;
        new_estimate = trapezoid_iterative_helper(start, end, steps, fp, old_estimate);
    }

    return new_estimate;
}

pub fn trapezoid_iterative_helper(
    start: f64,
    end: f64,
    steps: i64,
    fp: fn(f64) -> f64,
    old_estimate: f64,
) -> f64 {
    if steps == 1 {
        return (fp(start) + fp(end)) * (end - start) / 2.0;
    }

    let n = 1 << (steps - 2); // number of new points

    let h = (end - start) / (n as f64); // spacing of new points
    let x = start + h / 2.0; // coord of first new point
    let mut total = 0.0;

    for i in 0..n {
        total += fp(x + (i as f64) * h);
    }

    return (old_estimate + h * total) / 2.0;
}

pub fn romberg(start: f64, end: f64, _steps: i64, fp: fn(f64) -> f64) -> f64 {
    let max: usize = 21; // iterations

    let mut area = 0.0;
    let mut s = vec![0.0; max];

    for k in 1..max {
        let mut oldval = s[1];
        s[1] = trapezoid_iterative_helper(start, end, k as i64, fp, oldval);

        for i in 2..(k + 1) {
            let p = ((1 as i64) << (2 * (i - 1))) as f64;
            s[k] = (p * s[i - 1] - oldval) / (p - 1.0);
            oldval = s[i];
            s[i] = s[k];
        }

        if (area - s[k]).abs() < 1e-9 {
            return s[k];
        }

        area = s[k];
    }

    // max iterations reached; return what we have
    return s[max - 1];
}

pub fn simpsons(start: f64, end: f64, _steps: i64, fp: fn(f64) -> f64) -> f64 {
    let total = fp(start) + 4.0 * fp((start + end) / 2.0) + fp(end);
    return (end - start) / 6.0 * total;
}

pub fn simpsons38(start: f64, end: f64, _steps: i64, fp: fn(f64) -> f64) -> f64 {
    let total = fp(start)
        + 3.0 * fp((2.0 * start + end) / 3.0)
        + 3.0 * fp((2.0 * end + start) / 3.0)
        + fp(end);
    return (end - start) / 8.0 * total;
}

pub fn simpsons_composite(start: f64, end: f64, steps: i64, fp: fn(f64) -> f64) -> f64 {
    let steps = steps + (steps & 1);
    let h = (end - start) / (steps as f64);
    let mut totals = [0.0; 2];

    for i in 1..steps {
        totals[(i & 1) as usize] += fp(start + h * (i as f64));
    }

    let total = fp(start) + totals[0] * 2.0 + totals[1] * 4.0 + fp(end);

    return h * total / 3.0;
}

pub fn simpsons38_composite(start: f64, end: f64, steps: i64, fp: fn(f64) -> f64) -> f64 {
    let steps = steps * 3;

    let sections = steps / 3;

    let h = (end - start) / (steps as f64);
    let mut totals = [0.0; 4];

    for i in 0..sections {
        let mut x0 = start + 3.0 * (i as f64) * h;
        totals[0] += fp(x0);
        x0 += h;
        totals[1] += fp(x0);
        x0 += h;
        totals[2] += fp(x0);
        x0 += h;
        totals[3] += fp(x0);
    }

    let total = totals[0] + 3.0 * totals[1] + 3.0 * totals[2] + totals[3];

    return 3.0 * h / 8.0 * total;
}

pub fn bisect(fp: fn(f64) -> f64, mut x0: f64, mut x1: f64) -> f64 {
    let eps = 10e-8;
    let mut fx0 = fp(x0);
    let mut root = 0.0;

    while x1 - x0 > eps {
        root = (x0 + x1) / 2.0;
        let froot = fp(root);
        if froot * fx0 < 0.0 {
            x1 = root;
        } else {
            x0 = root;
            fx0 = froot;
        }
    }

    return root;
}

pub fn secant(fp: fn(f64) -> f64, mut x0: f64, mut x1: f64) -> f64 {
    let eps = 10e-8;

    let mut fx0 = fp(x0);
    while (x1 - x0).abs() > eps {
        let fx1 = fp(x1);
        let x2 = x1 - fx1 * ((x1 - x0) / (fx1 - fx0));
        x0 = x1;
        x1 = x2;
        fx0 = fx1;
    }

    return x1;
}

pub fn falsi(fp: fn(f64) -> f64, mut s: f64, mut t: f64) -> f64 {
    let mut side = 0;

    let mut fs = fp(s);
    let mut ft = fp(t);

    let m = 100; // max iterations

    let mut r = 0.0;
    let e = 5e-15; // eps

    for _ in 0..m {
        r = (fs * t - ft * s) / (fs - ft);
        if (t - s).abs() < (e * (t + s).abs()) {
            break;
        }

        let fr = fp(r);

        if fr * ft > 0.0 {
            t = r;
            ft = fr;
            if side == -1 {
                fs /= 2.0;
            }
            side = -1;
        } else if fs * fr > 0.0 {
            s = r;
            fs = fr;
            if side == 1 {
                ft /= 2.0;
            }
            side = 1;
        } else {
            break;
        }
    }
    return r;
}

pub fn ridder(fp: fn(f64) -> f64, mut a: f64, mut b: f64) -> f64 {
    let tol = 1.0e-9;

    let mut fa = fp(a);
    if fa == 0.0 {
        return a;
    }
    let mut fb = fp(b);
    if fb == 0.0 {
        return b;
    }

    // if fa*fb > 0.0: errror.err('Root is not bracketed')

    let mut xOld = 0.0;

    for i in 0..30 {
        // Compute the improved root x from Ridder's formula
        let c = 0.5 * (a + b);
        let fc = fp(c);
        let s = (fc * fc - fa * fb).sqrt();
        if s == 0.0 {
            return 0.0;
        }
        let mut dx = (c - a) * fc / s;
        if (fa - fb) < 0.0 {
            dx = -dx;
        }
        let x = c + dx;
        let fx = fp(x);
        // Test for convergence
        if i > 0 && (x - xOld).abs() < tol * x.abs().max(1.0) {
            return x;
        }
        xOld = x;
        // Re-bracket the root as tightly as possible
        if fc * fx > 0.0 {
            if fa * fx < 0.0 {
                b = x;
                fb = fx;
            } else {
                a = x;
                fa = fx;
            }
        } else {
            a = c;
            b = x;
            fa = fc;
            fb = fx;
        }
    }

    // too many iterations
    return 0.0;
}

#[cfg(test)]
mod tests {

    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn runall_integrate() {
        let algs = [
            trapezoid,
            trapezoid_iterative,
            romberg,
            simpsons_composite,
            simpsons38_composite,
            //  These two algorithms are not accurate enough
            // simpsons,
            // simpsons38
        ];

        for alg in algs.iter() {
            let r = alg(0.0, PI / 2.0, 100, f_integrate);
            let t = (r - 2.0).abs();
            assert!(t < 0.0001);
        }
    }

    fn f_integrate(x: f64) -> f64 {
        x.sin() + x.cos()
    }

    #[test]
    fn runall_roots() {
        let algs = [bisect, secant, falsi, ridder];

        for alg in algs.iter() {
            let r = alg(f_root, 1.0, 4.0);
            let t = (r - 2.0).abs();
            println!("t = {}, r = {}", t, r);
            assert!(t < 0.00001);
        }
    }

    fn f_root(x: f64) -> f64 {
        (x * x) - 4.0
    }
}
