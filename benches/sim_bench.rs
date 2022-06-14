use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rdiff_rs::{Simulation, SimulationParameters};

pub fn criterion_benchmark(c: &mut Criterion) {

    let mut params = SimulationParameters::coral_growth();
    let mut simulation = Simulation::new(100, 100, params);

    c.bench_function("sim_coral_100x100", |b| b.iter(|| (black_box(simulation.simulate()))));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);