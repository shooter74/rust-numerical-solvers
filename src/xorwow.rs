pub struct Xorwow {
    x: u32,
    y: u32,
    z: u32,
    w: u32,
    v: u32,
    d: u32,
}

impl Xorwow {
    pub fn new(seed: u32) -> Xorwow {
        let mut rng = Xorwow {
            x: 123456789,
            y: 362436069,
            z: 521288629,
            w: seed,
            v: 88675123,
            d: 6615241,
        };
        for _ in 0..10 {
            rng.next_u32();
        }
        rng
    }

    fn wrapping_add(&self, rhs: u32) -> u32 {
        let (sum, overflowed) = self.w.overflowing_add(rhs);
        if overflowed {
            let (sum, _) = sum.overflowing_add(1);
            sum
        } else {
            sum
        }
    }

    pub fn next_u32(&mut self) -> u32 {
        let t = self.x ^ (self.x >> 2);
        self.x = self.y;
        self.y = self.z;
        self.z = self.w;
        self.w = self.v;
        self.v = (self.v ^ (self.v << 4)) ^ (t ^ (t << 1));
        self.d += 362437;
        self.wrapping_add(self.v ^ self.d)
    }

    pub fn next_f64(&mut self) -> f64 {
        const NORM: f64 = 2.3283064365386963e-10;// 2^-32
        self.next_u32() as f64 * NORM
    }

    pub fn rand_vec(&mut self, n: usize) -> na::DVector<f64> {
        let mut vec = na::DVector::zeros(n);
        for i in 0..n {
            vec[i] = self.next_f64();
        }
        return vec;
    }
}
