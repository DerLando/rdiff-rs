use std::{mem::swap, path::Path};

use image::{ImageBuffer, Luma};



trait Empty {}
trait Initialized {}

#[derive(Debug)]
struct Cell {
    pub a: f64,
    pub b: f64,
    neighbor_indices: [usize; 8],
}

impl Default for Cell {
    fn default() -> Self {
        Self { a: 1.0, b: 0.0, neighbor_indices: Default::default() }
    }
}

impl Cell {
    pub fn neighbors<'a>(&'a self, grid: &'a Grid) -> [&Cell;8] {
        self.neighbor_indices.map(|i| &grid.cells[i])
    }
}

struct Grid {
    width: usize,
    height: usize,
    cells: Vec<Cell>
}

impl Grid {
    fn index_to_coords(index: usize, width: usize) -> (usize, usize) {
        let x = index % width;
        let y = index / width;

        (x, y)
    }

    fn coords_to_index(coords: (usize, usize), width: usize) -> usize {
        coords.1 * width + coords.0
    }

    fn neighbor_indices(index: usize, width: usize, height: usize) -> [usize; 8] {
        let pos = Grid::index_to_coords(index, width);
        let left = {
            if pos.0 == 0 {
                width - 1
            } else {
                pos.0 - 1
            }
        };
        let right = (pos.0 + 1) % width;
        let top = {
            if pos.1 == 0 {
                height - 1
            } else {
                pos.1 - 1
            }
        };
        let bot = (pos.1 + 1) % height;

        [
            Grid::coords_to_index((left, top), width),
            Grid::coords_to_index((pos.0, top), width),
            Grid::coords_to_index((right, top), width),
            Grid::coords_to_index((left, pos.1), width),
            Grid::coords_to_index((right, pos.1), width),
            Grid::coords_to_index((left, bot), width),
            Grid::coords_to_index((pos.0, bot), width),
            Grid::coords_to_index((right, bot), width),            
        ]

    }

    pub fn new(width: usize, height: usize) -> Self {
        let cells: Vec<Cell> = (0..(width * height)).map(|index| {
            let mut cell = Cell::default();
            cell.neighbor_indices = Grid::neighbor_indices(index, width, height);
            cell})
            .collect();
        Self { width, height, cells }
    }

    pub fn from_raw_parts(width: usize, height: usize, cells: Vec<Cell>) -> Self {
        Self {width, height, cells}
    }

    pub fn len(&self) -> usize {self.cells.len()}

    pub fn cell(&self, index: usize) -> &Cell {&self.cells[index]}

    pub fn cell_mut(&mut self, index: usize) -> &mut Cell {&mut self.cells[index]}

    pub fn cell_at_pos(&self, pos: (usize, usize)) -> &Cell {
        self.cell(Self::coords_to_index(pos, self.width))
    }

    pub fn cell_at_pos_mut(&mut self, pos: (usize, usize)) -> &mut Cell {
        self.cell_mut(Self::coords_to_index(pos, self.width))
    }
}

struct SimulationParameters {
    f: f64,
    k: f64,
    adj: f64,
    diag: f64,
    diff_a: f64,
    diff_b: f64,
    iter: usize
}

fn laplacian(cell: &Cell, grid: &Grid, params: &SimulationParameters) -> (f64, f64) {
    let neighbors = cell.neighbors(grid);
    let u = neighbors[1];
    let d = neighbors[6];
    let l = neighbors[3];
    let r = neighbors[4];
    
    let lu = neighbors[0];
    let ru = neighbors[2];
    let ld = neighbors[5];
    let rd = neighbors[7];
    
    (
        -cell.a + ((u.a + d.a + l.a + r.a) * params.adj) + ((lu.a + ru.a + ld.a + rd.a) * params.diag),
        -cell.b + ((u.b + d.b + l.b + r.b) * params.adj) + ((lu.b + ru.b + ld.b + rd.b) * params.diag)
    )
}
struct Simulation {
    current_grid: Grid,
    next_grid: Grid,
    iter: usize,
    params: SimulationParameters
}

impl Simulation {
    pub fn new(width: usize, height: usize, params: SimulationParameters) -> Self {
        Self { current_grid: Grid::new(width, height), next_grid: Grid::new(width, height), iter: 0, params }
    }

    pub fn seed_cell(&mut self, pos: (usize, usize)) {
        self.current_grid.cell_at_pos_mut(pos).b = 1.0;
    }

    pub fn seed_dot(&mut self, pos: (usize, usize)) {
        let positions = [
            (pos.0, pos.1 - 1),
            (pos.0 - 1, pos.1),
            pos,
            (pos.0 + 1, pos.1),
            (pos.0, pos.1 + 1)
        ];

        for position in positions {
            self.seed_cell(position);
        }
    }

    fn advance(&mut self) {
        for i in 0..self.current_grid.len() {
            let cell = self.current_grid.cell(i);
            let (lap_a, lap_b) = laplacian(cell, &self.current_grid, &self.params);
            let abb = cell.a * cell.b * cell.b;

            let a = cell.a + (self.params.diff_a * lap_a - abb + self.params.f * (1.0 - cell.a));
            let b = cell.b + (self.params.diff_b * lap_b + abb - (self.params.k + self.params.f) * cell.b );

            self.next_grid.cell_mut(i).a = a;
            self.next_grid.cell_mut(i).b = b;
        }

        // swap next and current
        swap(&mut self.current_grid, &mut self.next_grid);

        // raise iteration count
        self.iter += 1;
    }

    pub fn simulate(&mut self) {
        for _ in 0..self.params.iter {
            self.advance();
        }
    }

    pub fn save_to_image(&self, filename: &impl AsRef<Path>) {
        fn luma(cell: &Cell) -> Luma<u8>  {
            Luma([(cell.a * 255.0) as u8])
        }
        let mut img = ImageBuffer::from_fn(self.current_grid.width as u32, self.current_grid.height as u32, |x, y| {luma(self.next_grid.cell_at_pos((x as usize, y as usize)))});

        img.save(filename);
    }
}


#[cfg(test)]
mod tests {
    use crate::{Grid, Simulation};

    #[test]
    fn it_works() {
        let grid = Grid::new(1000, 300);

        assert_eq!(Grid::coords_to_index((0, 0), 1000), 0);
        assert_eq!(Grid::coords_to_index((0, 1), 1000), 1000);
        assert_eq!(Grid::coords_to_index((0, 2), 1000), 2000);
        assert_eq!(Grid::coords_to_index((999, 2), 1000), 2999);

    }

    #[test]
    fn should_produce_image() {
        let params = crate::SimulationParameters { f: 0.0545, k: 0.062, adj: 0.2, diag: 0.05, diff_a: 1.0, diff_b: 0.5, iter: 2000 };
        let mut simulation = Simulation::new(100, 100, params);

        simulation.seed_dot((25, 25));
        // simulation.seed_dot((75, 25));
        // simulation.seed_dot((75, 75));
        // simulation.seed_dot((25, 75));
        // simulation.seed_dot((50, 50));


        simulation.simulate();

        simulation.save_to_image(&"test.png");
    }
}