#![feature(test)]

extern crate test;

pub mod functions {

    use std::f64;

    #[inline]
    pub fn ellips_func(x: &[f64]) -> f64 {
        let mut f: f64 = 0.0;
        let mut n: f64 = 0.0;
        let med: f64 = x.len() as f64 - 1.0;
        for i in 0..x.len() {
            f += f64::powf(1000000.0, n / med) * x[i] * x[i];
            n += 1.0;
        }
        f
    }

    pub fn bent_cigar_func(x: &[f64]) -> f64 {
        let mut f: f64 = 0.0;
        let f_0: f64 = x[0] * x[0];
        for i in 1..x.len() {
            f += x[i] * x[i];
        }
        f_0 + 1000000.0 * f
    }

    pub fn discuss_func(x: &[f64]) -> f64 {
        let mut f: f64 = 0.0;
        let f_0: f64 = x[0] * x[0];
        for i in 1..x.len() {
            f += x[i] * x[i];
        }
        1000000.0 * f_0 + f
    }

    pub fn rosenbrock_func(x: &[f64]) -> f64 {
        let mut f: f64 = 0.0;
        for i in 0..x.len() - 1 {
            f += 100.0 * (x[i] * x[i] - x[i + 1]) * (x[i] * x[i] - x[i + 1])
                + (x[i] - 1.0) * (x[i] - 1.0);
        }
        f
    }

    pub fn ackley_func(x: &[f64]) -> f64 {
        let mut sum1: f64 = 0.0;
        let mut sum2: f64 = 0.0;
        let dimension: f64 = x.len() as f64;
        for i in 0..x.len() {
            sum1 += x[i] * x[i];
            sum2 += f64::cos(2.0 * f64::consts::PI * x[i]);
        }
        -20.0 * f64::exp(-0.2 * (sum1 / dimension).sqrt()) - f64::exp(sum2 / dimension)
            + 20.0
            + f64::exp(1.0)
    }

    pub fn weierstrass_func(x: &[f64]) -> f64 {
        let kmax = 20;
        let a = 0.5;
        let b = 3.0;

        let mut sum1: f64 = 0.0;

        for i in 0..x.len() {
            let mut k: f64 = 0.0;
            for _ in 0..kmax + 1 {
                sum1 += f64::powf(a, k)
                    * f64::cos(2.0 * f64::consts::PI * f64::powf(b, k) * (x[i] + 0.5));
                k += 1.0;
            }
        }

        let mut sum2: f64 = 0.0;
        let mut k: f64 = 0.0;
        for _ in 0..kmax + 1 {
            sum2 += f64::powf(a, k) * f64::cos(2.0 * f64::consts::PI * f64::powf(b, k) * 0.5);
            k += 1.0;
        }
        sum1 - x.len() as f64 * sum2
    }

    pub fn griewank_func(x: &[f64]) -> f64 {
        let mut sum1: f64 = 0.0;
        let mut sum2: f64 = 1.0;
        let mut _i: f64 = 1.0;
        for i in 0..x.len() {
            sum1 += x[i] * x[i];
            sum2 *= f64::cos(x[i] / _i.sqrt());
            _i += 1.0;
        }
        1.0 + sum1 / 4000.0 - sum2
    }

    pub fn rastrigin_func(x: &[f64]) -> f64 {
        let mut f: f64 = 0.0;
        for i in 0..x.len() {
            f += x[i] * x[i] - 10.0 * f64::cos(2.0 * f64::consts::PI * x[i]) + 10.0;
        }
        f
    }

    pub fn schwefel_func(x: &[f64]) -> f64 {
        let mut sum: f64 = 0.0;
        for i in 0..x.len() {
            sum += x[i] * f64::sin(f64::abs(x[i]).sqrt())
        }
        418.9829 * x.len() as f64 - sum
    }
}

#[cfg(test)]
mod tests_functions {
    use functions;
    use test::Bencher;

    #[bench]
    fn bench_bent_cigar_func(b: &mut Bencher) {
        let _array: [f64; 30] = [0.0; 30];
        b.iter(|| functions::bent_cigar_func(&_array[0..30]));
    }

    #[bench]
    fn bench_ackley_func(b: &mut Bencher) {
        let _array: [f64; 30] = [0.0; 30];
        b.iter(|| functions::ackley_func(&_array[0..30]));
    }

    #[bench]
    fn bench_weierstrass_func(b: &mut Bencher) {
        let _array: [f64; 30] = [0.0; 30];
        b.iter(|| functions::weierstrass_func(&_array[0..30]));
    }

    #[bench]
    fn bench_griewank_func(b: &mut Bencher) {
        let _array: [f64; 30] = [0.0; 30];
        b.iter(|| functions::griewank_func(&_array[0..30]));
    }

    #[bench]
    fn bench_rastrigin_func(b: &mut Bencher) {
        let _array: [f64; 30] = [0.0; 30];
        b.iter(|| functions::rastrigin_func(&_array[0..30]));
    }

    #[bench]
    fn bench_schwefel_func(b: &mut Bencher) {
        let _array: [f64; 30] = [0.0; 30];
        b.iter(|| functions::schwefel_func(&_array[0..30]));
    }

}
