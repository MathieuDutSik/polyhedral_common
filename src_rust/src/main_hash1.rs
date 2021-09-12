pub struct ShiftXorHasher {
    state: u64,
}

impl std::hash::Hasher for ShiftXorHasher {
    fn write(&mut self, bytes: &[u8]) {
        for &byte in bytes {
            self.state = self.state.rotate_left(8) ^ u64::from(byte);
        }
    }
    
    fn finish(&self) -> u64 {
        eprintln!("hasher called: {:#010x}", self.state);
        self.state
    }
}

pub struct BuildShiftXorHasher;

impl std::hash::BuildHasher for BuildShiftXorHasher {
    type Hasher = ShiftXorHasher;
    fn build_hasher(&self) -> ShiftXorHasher {
        ShiftXorHasher { state: 0 }
    }
}

fn main() {
    let mut hm = std::collections::HashMap::with_hasher(BuildShiftXorHasher);
    hm.insert(3, 4);
    hm.insert(5, 6);
    hm.insert(7, 8);
}
