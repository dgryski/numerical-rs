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

#[cfg(test)]
mod tests {

    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn runall() {
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
            let r = alg(0.0, PI/2.0, 100, f);
            let t = (r - 2.0).abs();
            assert!(t < 0.0001);
        }
    }

    fn f(x:f64) -> f64 { x.sin() + x.cos() }
}
