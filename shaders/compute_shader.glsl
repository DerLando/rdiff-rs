#version 430

layout(local_size_x = 1) in;

uniform vec2 texture_size;

struct Cell {
    float a;
    float b;
}

// Input buffer
layout(std430, binding=0) buffer cells_in {
    Cell cells[];
} In;

// Output buffer
layout(std430, binding=1) buffer cells_out {
    Cell cells[];
} Out;

void main() {
    int curCellIndex = int(gl_GlobalInvocationID);

    Cell in_cell = In.cells[curCellIndex];
}