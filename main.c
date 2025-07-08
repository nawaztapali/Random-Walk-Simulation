/*
 * Assessment 3 for PHY2027
 * Author: Mohammed Nawaz Tapali
 * Date: 12/12/2024
 */

/* Program description: This program does 1D, 2D, and 3D lattice simulations of random and self-avoiding random walks,
 * modeling particle motion on a lattice. Random walks model particle diffusion in which the direction of each step is a
 * random process, determined by the product of the diffusion constant, step size, and time interval. Self-avoiding random
 * walks impose some constraint: the particle is forbidden from returning to previously occupied positions. It is especially
 * suitable for physical systems, such as chains of polymers. The program calculates mean squared displacement and mean absolute
 * displacement, storing them in text files that can be analyzed. It also exports the paths of the particle in 1D, 2D, and 3D format,
 * which enables the visualization of the simulated motion. By investigating these simulations, the program gives insight into random
 * processes, diffusion, and constrained movement in complex systems.
 */




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

#define NUM_STEPS 10   // Number of steps for the random walk
#define NUM_RUNS 10  // Number of independent runs for the simulation
#define GRID_SIZE (25) // Size of the lattice
#define OFFSET (GRID_SIZE / 2)  // Offset to center the grid

// Function to generate random numbers
double generate_random_points() {
    return (double)rand() / RAND_MAX;
}

// Function to perform random walk in 1 dimension
// Parameters:
// steps: Number of steps to simulate
// dx: Step size (distance moved per step)
// dt: Time step (duration of each step)
// D: Diffusion constant
// x: Array to store the x-coordinate of the particle at each step
void random_walk_1D (int steps, float dx, float dt, float D, double *x) {
    double p = D * dt / pow(dx,2); //Probability


    for (int i = 1; i < NUM_STEPS; i++) {
        double r = generate_random_points();
        if (r < p) {
            x[i] = x[i - 1] + dx;    // movement in positive x-direction
        } else if (r < 2 * p){
            x[i] = x[i - 1] - dx;    // movement in negative x-direction
        } else {
            x[i] = x[i - 1];         // No movement
        }
    }
}



// Function to perform random walk in 2 dimension
// Parameters:
// steps: Number of steps to simulate
// dx: Step size (distance moved per step)
// dt: Time step (duration of each step)
// D: Diffusion constant
// x, y: Arrays to store the x and y coordinates of the particle at each step
void random_walk_2D (int steps, float dx, float dt, float D, double *x, double *y ) {
    double p = D * dt / pow(dx, 2);

    for (int i = 1; i < NUM_STEPS; i++) {
        double rx = generate_random_points();
        double ry = generate_random_points();
        if (rx < p) {
            x[i] = x[i - 1] + dx;       // movement in positive x-direction
        } else if (rx < 2 * p) {
            x[i] = x[i - 1] - dx;       // movement in negative x-direction
        } else {
            x[i] = x[i - 1];            // No movement
        }

        if (ry < p){
            y[i] = y[i - 1] + dx;       // movement in positive y-direction
        } else if (ry < 2 * p) {
            y[i] = y[i - 1] - dx;       // movement in negative y-direction
        } else {
            y[i] = y[i - 1];            // no movement
        }
    }

}



// Function to perform random walk in 3 dimension
// Parameters:
// steps: Number of steps to simulate
// dx: Step size (distance moved per step)
// dt: Time step (duration of each step)
// D: Diffusion constant
// x, y, z: Arrays to store the x, y, and z coordinates of the particle at each step
void random_walk_3D(int steps, float dx, float dt, float D, double *x, double *y, double *z) {
    double p = D * dt / pow(dx, 2);

    for (int i = 1; i < NUM_STEPS; i++) {
        double rx = generate_random_points();
        double ry = generate_random_points();
        double rz = generate_random_points();
        if (rx < p) {
            x[i] = x[i - 1] + dx;          // movement in positive x-direction
        } else if (rx < 2 * p) {
            x[i] = x[i - 1] - dx;          // movement in negative x-direction
        } else {
            x[i] = x[i - 1];               // no movement
        }

        if (ry < p) {
            y[i] = y[i - 1] + dx;          // movement in positive y-direction
        } else if (ry < 2 * p) {
            y[i] = y[i - 1] - dx;          // movement in negative y-direction
        } else {
            y[i] = y[i - 1];               // no movement
        }

        if (rz < p) {
            z[i] = z[i - 1] + dx;          // movement in positive z-direction
        } else if (rz < 2 * p){
            z[i] = z[i - 1] - dx;          // movement in negative z-direction
        } else {
            z[i] = z[i -1];                // no movement
        }
    }

}



// Function to check if a move is valid or not in 2 dimension
// Parameters:
// x, y: Proposed new coordinates
// grid: 2D array representing the grid to check boundaries
// Returns: 1 if the move is valid, 0 otherwise
int is_valid_move_2D(int x, int y, int grid[GRID_SIZE][GRID_SIZE]) {
    if (x < 0 || x >= GRID_SIZE || y < 0 || y >= GRID_SIZE) { // Checks if the particle is beyond the edge. If it is not there, it returns 1
        return 0;
    } else {
        return 1;
    }
}



// Function to check if a move is valid or not in 3 dimension
// Parameters:
// x, y, z: Proposed new coordinates
// visited: 3D array tracking visited cells
// Returns: 1 if the move is valid and the cell is unvisited, 0 otherwise
int is_valid_move_3D(int x, int y, int z, int visited[GRID_SIZE][GRID_SIZE][GRID_SIZE]) {
    if (x < 0 || x >= GRID_SIZE || y < 0 || y >= GRID_SIZE || z < 0 || z >= GRID_SIZE) {
        return 0; // it gives zero when the particle is out of bounds (invalid move if out of bounds), 0 otherwise
    }
    return !visited[x][y][z]; // Valid if not visited
}


// Function to perform self avoiding in 2D

void self_avoiding_random_walk_2D (int steps, float dx, float dt, float D, int grid_2D[GRID_SIZE][GRID_SIZE]) {
    // Array to keep track of whether each grid cell has been visited
    int visited[GRID_SIZE][GRID_SIZE] = {0}; //initialize an array to track the visited point; initialize all cells as not visited (0).


    // Start the particle at the center of the grid
    int grid_x = (int)OFFSET, grid_y = (int)OFFSET; // Start in the middle of the grid, so setting/initializing the movement from center
    grid_2D[grid_x ][grid_y] = 1;  // Marking the starting point which is 1 in this case. It could be any other number
    visited[grid_x ][grid_y] = 1; // Mark starting point as visited


    // Loop through the number of steps specified
    for (int i = 1; i < steps; i++) {
        int move_x[4] = {1, -1, 0, 0}; //0 and 0 represents restricted movement in top and bottom due to movement only in x
        int move_y[4] = {0, 0, 1, -1}; //Here, also 0 and 0 represents restricted movement in right and left due to movement only in y

        // Arrays to track which moves are valid
        int valid_moves[4] = {0};  // Initialize all moves as invalid (0).
        int valid_count = 0;       // Counter for the valid moves

        // Looping to check all possible moves
        for (int j = 0; j < 4; j++) {
            int new_x = grid_x + move_x[j];   // Updates the value of x by the allowed movements int move_x array
            int new_y = grid_y + move_y[j];   // Updates the value of y by the allowed movements int move_y array

            // Check if the move is valid(within bounds and not already visited)
            int result = is_valid_move_2D(new_x, new_y, visited);


            if (result == 1) {      // if the move is valid, it marks this move as valid and increment the counts
                valid_moves[j] = 1; // this tells what are the valid moves from move_x and move_y arrays
                valid_count++;
            }
        }


        // If no valid moves, terminate early
        if (valid_count == 0) {
            printf("Walk terminated early at step %d due to no valid moves.\n", i);
            break; // Exit the loop
        }

        // Randomly selects a valid move
        int selected_move = -1;  // Initialize the selected move as invalid
        while (selected_move == -1) { // keep choosing until a valid move is selected.
            double r = generate_random_points(); // Generates a random number between 0 and 1
            int candidate_move = (int)(r * 4);   //This is the suggested move between 0 and less than 4 by scaling them to one of the possible move

            // calculate the candidate move's new position
            int c_x = grid_x + move_x[candidate_move];
            int c_y = grid_y + move_y[candidate_move];

            // If the move is valid and not already visited, it selects the move
            if (valid_moves[candidate_move] == 1 && visited[c_x][c_y] == 0) {
                selected_move = candidate_move; // Assign the chosen move
            }
        }

        int previous_grid_value = grid_2D[grid_x][grid_y];

        // Update the partilce's position based on the selected move
        grid_x += move_x[selected_move]; // Updates the x-coordinate
        grid_y += move_y[selected_move]; //updates the y-coordinate

        // Mark the new position as visited
        visited[grid_x][grid_y] = 1; // update the visited array
        grid_2D[grid_x][grid_y] = previous_grid_value + 1;    // Mark the grid to visualize the particle's path
    }

}




// 3D Self-Avoiding Random Walk Function
void self_avoiding_random_walk_3D(int steps, float dx, float dt, float D, int grid[GRID_SIZE][GRID_SIZE][GRID_SIZE]) {

    int visited[GRID_SIZE][GRID_SIZE][GRID_SIZE] = {0}; // Initialize a 3D array to track visited points

    int grid_x = (int)OFFSET, grid_y = (int)OFFSET, grid_z = (int)OFFSET; // Start in the center
    grid[grid_x][grid_y][grid_z] = 5;  // Initialize the starting point
    visited[grid_x][grid_y][grid_z] = 1; // Mark starting point as visited

    // Begin the simulation loop for the given number of steps
    for (int i = 1; i < steps; i++) {

        int move_x[6] = {1, -1, 0, 0, 0, 0}; // Movement in x-direction
        int move_y[6] = {0, 0, 1, -1, 0, 0}; // Movement in y-direction
        int move_z[6] = {0, 0, 0, 0, 1, -1}; // Movement in z-direction

        // Array to track valid moves and a counter for how many valid moves are available
        int valid_moves[6] = {0};
        int valid_count = 0;

        // Check all possible moves to ensure they don't revisit a point or go out of bounds
        for (int j = 0; j < 6; j++) {
            int new_x = grid_x + move_x[j];
            int new_y = grid_y + move_y[j];
            int new_z = grid_z + move_z[j];

            if (is_valid_move_3D(new_x, new_y, new_z, visited)) {
                valid_moves[j] = 1; // Mark move as valid
                valid_count++;
            }
        }

        // If no valid moves, terminate early
        if (valid_count == 0) {
            printf("Walk terminated early at step %d due to no valid moves.\n", i);
            break;
        }

        int selected_move = -1;
        while (selected_move == -1) {
            double r = generate_random_points();
            int candidate_move = (int)(r * 6); // Randomly choose a move among the 6 possible directions
            int c_x = grid_x + move_x[candidate_move];
            int c_y = grid_y + move_y[candidate_move];
            int c_z = grid_z + move_z[candidate_move];
            if (valid_moves[candidate_move] == 1 && visited[c_x][c_y][c_z] == 0) {
                selected_move = candidate_move;
            }
        }

        // Save the previous value at the grid position for incrementing
       int last_grid_value = grid[grid_x][grid_y][grid_z];

        // Update position to the new coordinates
        grid_x += move_x[selected_move];
        grid_y += move_y[selected_move];
        grid_z += move_z[selected_move];

        // Mark the new position as visited
        visited[grid_x][grid_y][grid_z] = 1;

        // Update the grid to indicate the new position
        grid[grid_x][grid_y][grid_z] = last_grid_value + 1;    // Update the grid to indicate the new position
    }
}




// Export results to a file. This function exports the results of the walk in 1D, 2D, or 3D format
void write_to_file(const char *filename, double *x, double *y, double *z, int steps) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        printf("Error: Could not open file %s for writing.\n", filename);
        return;
    }

    for (int i = 0; i < steps; i++) {
        if (z) {
            fprintf(file, "%f %f %f\n", x[i], y[i], z[i]);  // 3D data
        } else if (y) {
            fprintf(file, "%f %f\n", x[i], y[i]);          // 2D data
        } else {
            fprintf(file, "%f\n", x[i]);                  // 1D data
        }
    }

    fclose(file);
}



// Export results to a file. This function exports the results of a 2D grid to a file
void write_to_file_2D(const char *filename, int grid[GRID_SIZE][GRID_SIZE], int steps) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        printf("Error: Could not open file %s for writing.\n", filename);
        return;
    }

    // Iterate through each step of the walk
    int current_step = 1;
    while(current_step <= steps) {
        for (int i = 0; i < GRID_SIZE; i++) {
            for (int j = 0; j < GRID_SIZE; j++) {
                if (grid[i][j] == current_step) {   // Check if the cell corresponds to the current step
                    fprintf(file, "%d %d\n", i, j); // Write the coordinates

                }
            }
        }
        current_step++; // Move to the next step
    }

    fclose(file);
}


// Function to export results of a 3D grid to a fil
void write_to_file_3D(const char *filename, int grid[GRID_SIZE][GRID_SIZE][GRID_SIZE], int steps) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        printf("Error: Could not open file %s for writing.\n", filename);
        return;
    }

    // Iterate through each step of the walk
    int current_step = 1;
    while(current_step <= steps) {
        for (int i = 0; i < GRID_SIZE; i++) {
            for (int j = 0; j < GRID_SIZE; j++) {
                for (int k = 0; k < GRID_SIZE; k++) {
                    if (grid[i][j][k] == current_step) {       // Check if the cell corresponds to the current ste
                        fprintf(file, "%d %d %d\n", i, j, k);
                    }
                }
            }
        }
        current_step++;
    }

    fclose(file);
}









// Function to write results of a 2D self-avoiding walk simulation to a file
void write_results_to_file_2D(int num_steps, int num_runs, const char *filename) {
    // Open the specified file for writing
    FILE *file = fopen(filename, "w");

    const float dt = 0.1;  // Define the time step for simulation data

    // Loop through each simulation step
    for (int step = 1; step <= num_steps; step++) {
        double mean_r = 0.0, mean_r_squared = 0.0; // Initialize variables to calculate mean results

        // Perform multiple runs for the same step to calculate averages
        for (int run = 0; run < num_runs; run++) {
            // Initialize a 2D grid to track visited points during the walk
            bool visited[GRID_SIZE][GRID_SIZE] = {false};
            int x = OFFSET, y = OFFSET; // Start at the central point of the grid

            visited[x][y] = true;  // Mark the starting position as visited

            // Perform a single self-avoiding walk up to the current number of steps
            for (int t = 1; t <= step; t++) {
                // Define possible moves in the 2D grid
                int possible_moves[4][2] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
                int valid_moves[4], valid_count = 0; // Track valid moves and their count

                // Determine which moves are valid
                for (int i = 0; i < 4; i++) {
                    int nx = x + possible_moves[i][0]; // Calculate new x-coordinate
                    int ny = y + possible_moves[i][1]; // Calculate new y-coordinate

                    // Check if the new position is within bounds and unvisited
                    if (nx >= 0 && nx < GRID_SIZE && ny >= 0 && ny < GRID_SIZE && !visited[nx][ny]) {
                        valid_moves[valid_count++] = i; // Add valid move index to the list
                    }
                }

                // Terminate the walk if no valid moves are available
                if (valid_count == 0) {
                    break;
                }

                // Select a random valid move from the list
                int selected_move = valid_moves[rand() % valid_count];
                x += possible_moves[selected_move][0]; // Update x-coordinate
                y += possible_moves[selected_move][1]; // Update y-coordinate
                visited[x][y] = true;  // Mark the new position as visited
            }

            // Calculate displacement from the starting position
            int displacement_x = x - OFFSET, displacement_y = y - OFFSET;
            double r_squared = pow(displacement_x, 2) + pow(displacement_y, 2); // Squared displacement
            double r_abs = sqrt(r_squared); // Absolute displacement

            // Accumulate results for this run
            mean_r_squared += r_squared;
            mean_r += r_abs;
        }

        // Calculate averages for this step over all runs
        mean_r_squared /= num_runs;
        mean_r /= num_runs;

        // Write the results for this step to the file
        fprintf(file, "%f %f %f\n", step * dt, mean_r, mean_r_squared);
    }

    fclose(file); // Close the file after writing
    printf("Results saved to %s\n", filename); // Notify user of successful save
}




// Function to write results of a 3D self-avoiding walk simulation to a file
void write_results_to_file_3D(int num_steps, int num_runs, const char *filename) {
    // Open the specified file for writing
    FILE *file = fopen(filename, "w");

    const float dt = 0.1;  // Define the time step for simulation data

    // Loop through each simulation step
    for (int step = 1; step <= num_steps; step++) {
        double mean_r = 0.0, mean_r_squared = 0.0; // Initialize variables to calculate mean results

        // Perform multiple runs for the same step to calculate averages
        for (int run = 0; run < num_runs; run++) {
            // Initialize a 3D grid to track visited points during the walk
            bool visited[GRID_SIZE][GRID_SIZE][GRID_SIZE] = {false};
            int x = OFFSET, y = OFFSET, z = OFFSET; // Start at the central point of the grid

            visited[x][y][z] = true;  // Mark the starting position as visited

            // Perform a single self-avoiding walk up to the current number of steps
            for (int t = 1; t <= step; t++) {
                // Define possible moves in the 3D grid
                int possible_moves[6][3] = {
                    {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
                int valid_moves[6], valid_count = 0; // Track valid moves and their count

                // Determine which moves are valid
                for (int i = 0; i < 6; i++) {
                    int nx = x + possible_moves[i][0]; // Calculate new x-coordinate
                    int ny = y + possible_moves[i][1]; // Calculate new y-coordinate
                    int nz = z + possible_moves[i][2]; // Calculate new z-coordinate

                    // Check if the new position is within bounds and unvisited
                    if (nx >= 0 && nx < GRID_SIZE &&
                        ny >= 0 && ny < GRID_SIZE &&
                        nz >= 0 && nz < GRID_SIZE &&
                        !visited[nx][ny][nz]) {
                        valid_moves[valid_count++] = i; // Add valid move index to the list
                    }
                }

                // Terminate the walk if no valid moves are available
                if (valid_count == 0) {
                    break;
                }

                // Select a random valid move from the list
                int selected_move = valid_moves[rand() % valid_count];
                x += possible_moves[selected_move][0]; // Update x-coordinate
                y += possible_moves[selected_move][1]; // Update y-coordinate
                z += possible_moves[selected_move][2]; // Update z-coordinate
                visited[x][y][z] = true;  // Mark the new position as visited
            }

            // Calculate displacement from the starting position
            int displacement_x = x - OFFSET;
            int displacement_y = y - OFFSET;
            int displacement_z = z - OFFSET;

            double r_squared = pow(displacement_x, 2) + pow(displacement_y, 2) + pow(displacement_z, 2); // Squared displacement
            double r_abs = sqrt(r_squared); // Absolute displacement

            // Accumulate results for this run
            mean_r_squared += r_squared;
            mean_r += r_abs;
        }

        // Calculate averages for this step over all runs
        mean_r_squared /= num_runs;
        mean_r /= num_runs;

        // Write the results for this step to the file
        fprintf(file, "%f %f %f\n", step * dt, mean_r, mean_r_squared);
    }

    fclose(file); // Close the file after writing
    printf("Results saved to %s\n", filename); // Notify user of successful save
}





int main(void) {

    const float DX = 1;   //step size [m]
    const float DT = 0.1; //time step [s]
    const float D = 2;    //diffusion constant [m^2/s]
    double x_coordinate[NUM_STEPS] = {0};
    double y_coordinate[NUM_STEPS] = {0};
    double z_coordinate[NUM_STEPS] = {0};

    //initialization of lattice points in 2D
    int grid_2D[GRID_SIZE][GRID_SIZE];
    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            grid_2D[i][j] = 0;
        }
    }

    // Initialize the 3D lattice
    int grid_3D[GRID_SIZE][GRID_SIZE][GRID_SIZE];
    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            for (int k = 0; k < GRID_SIZE; k++) {
                grid_3D[i][j][k] = 0;
            }
        }
    }

    // Seed the random number generator
    srand(time(NULL));
    rand();


    /* This section of the code is calling all the functions that execute the simulation of random walks
     * in different dimensions by calling the prototype functions which have created earlier in the code
     * before the main function. The functions which are called here are random_walk_1D, random_walk_2D,
     * random_walk_3D, self_avoiding_random_walk_2D, and self_avoiding_random_walk_3D.
     */

    // Performing 1D random walk
    random_walk_1D(NUM_STEPS, DX, DT, D, x_coordinate);
    // Performing 2D random walk
    random_walk_2D(NUM_STEPS, DX, DT, D, x_coordinate, y_coordinate);
    // Performing 3D random walk
    random_walk_3D(NUM_STEPS, DX, DT, D, x_coordinate, y_coordinate, z_coordinate);
    // Performing 2D self-avoiding random walk
    self_avoiding_random_walk_2D(NUM_STEPS, DX, DT, D, grid_2D);
    // Performing 2D the self-avoiding random walk
    self_avoiding_random_walk_3D(NUM_STEPS, DX, DT, D, grid_3D);


    /* This section of the codes provides a way to visualize the simulation of random walks. This is done
     * by calling the functions which are devised to write the results. The first functions names is write_to_file.
     * prints, and it prints the coordinates in to the file. For 1D, it prints x-coordinates. For 2D, it prints
     * x-coordinates and y-coordinates. For 3D,it prints x-coordinates, y-coordinates and z-coordinates. The other
     * two functions whose names are write_to_file_2D and write_to_file_3D print the x and y coordinates and x, y, and
     * z coordinates respectively form the center because they start form the the center of the lattice.
     */

    //Writing the x-coordinates to the file
    write_to_file("random_walk_1D.txt", x_coordinate, NULL, NULL, NUM_STEPS);

    //Writing the x-coordinates and y-coordinates to the file
    write_to_file("random_walk_2D.txt", x_coordinate, y_coordinate, NULL, NUM_STEPS);

    //Writing the x-coordinates, y-coordinates and z-coordinates to the file
    write_to_file("random_walk_3D.txt", x_coordinate, y_coordinate, z_coordinate, NUM_STEPS);

    //Writing the x-coordinates and y-coordinates to the file from the center
    write_to_file_2D("self_avoiding_random_walk_2D.txt", grid_2D, NUM_STEPS);

    //Writing the x-coordinates,y-coordinates and z-coordinates to the file from the center
    write_to_file_3D("self_avoiding_random_walk_3D.txt", grid_3D, NUM_STEPS);



    /* This section of the code particularly deals with the calculation of average of position <x> and
     * the average of position squared <x�>. This is dealt separately with each dimensions. Firstly, arrays
     * has been created to store the calculation and then they are initialized. Then, multiple random walks
     * are executed. Simultaneously, the average of x and average of x-squared are calculated. Then, it is
     * printed in the file for the plotting purposes.
     */


    // Accumulators for averages for 1D random walks
    double x_squared_avg[NUM_STEPS], x_abs_avg[NUM_STEPS], x[NUM_STEPS];
    const float dt = 0.1;  // Time step
    const char *filename = "1D_results.txt";  // File to save results

    // Initialize averages to zero
    for (int i = 0; i < NUM_STEPS; i++) {
        x_squared_avg[i] = 0;
        x_abs_avg[i] = 0;
    }

    // Loop over multiple random walks
    for (int run = 0; run < NUM_RUNS; run++) {
        x[0] = 0;
        // Perform a single random walk
        random_walk_1D(NUM_STEPS, DX, DT, D, x);

            // Accumulate results for each time step
            for (int t = 0; t < NUM_STEPS; t++) {
                x_squared_avg[t] += pow(x[t], 2);  // Sum of squared displacements
                x_abs_avg[t] += fabs(x[t]);       // Sum of absolute displacements
        }
    }

    // Compute averages
    for (int t = 0; t < NUM_STEPS; t++) {
        x_squared_avg[t] /= NUM_RUNS;
        x_abs_avg[t] /= NUM_RUNS;
    }

    // Save results to a file
    FILE *file = fopen(filename, "w");


    // Write data to file: time (T), <x�>, <|x|>
    for (int t = 1; t < NUM_STEPS; t++) {  // Start from t = 1
        fprintf(file, "%f %f %f\n", t * dt, x_squared_avg[t], x_abs_avg[t]);
    }

    fclose(file);
    printf("Results saved to %s\n", filename);





    // Accumulators for averages for 2D random walks
    double r_squared_avg[NUM_STEPS], r_abs_avg[NUM_STEPS];
    double x_2D[NUM_STEPS], y_2D[NUM_STEPS];  // Ensure 'x' is not redeclared if already defined earlier
    const char *filename_2D = "results_2D.txt";  // Correct declaration

    // Initialize averages to zero
    for (int i = 0; i < NUM_STEPS; i++) {
        r_squared_avg[i] = 0;
        r_abs_avg[i] = 0;
    }

    // Loop over multiple random walks
    for (int run = 0; run < NUM_RUNS; run++) {
        x_2D[0] = 0;
        y_2D[0] = 0;

        // Perform a single 2D random walk
        random_walk_2D(NUM_STEPS, DX, DT, D, x_2D, y_2D);

        // Accumulate results for each time step
        for (int t = 0; t < NUM_STEPS; t++) {
            double r_squared = pow(x_2D[t], 2) + pow(y_2D[t], 2);
            double r_abs = sqrt(r_squared);
            r_squared_avg[t] += r_squared;
            r_abs_avg[t] += r_abs;
        }
}

    // Compute averages
    for (int t = 0; t < NUM_STEPS; t++) {
        r_squared_avg[t] /= NUM_RUNS;
        r_abs_avg[t] /= NUM_RUNS;
    }

    // Save results to a file
    FILE *file_2D = fopen(filename_2D, "w");


    // Write data to file: time (T), <r�>, <|r|>
    for (int t = 1; t < NUM_STEPS; t++) {  // Start from t = 1
        fprintf(file_2D, "%f %f %f\n", t * dt, r_squared_avg[t], r_abs_avg[t]);
    }

    fclose(file_2D);
    printf("Results saved to %s\n", filename_2D);





    // Accumulators for averages for 3D random walks
    double r_squared_avg_3D[NUM_STEPS], r_abs_avg_3D[NUM_STEPS];
    double x_3D[NUM_STEPS], y_3D[NUM_STEPS], z_3D[NUM_STEPS];
    const char *filename_3D = "results_3D.txt";

    // Initialize averages to zero
    for (int i = 0; i < NUM_STEPS; i++) {
        r_squared_avg_3D[i] = 0;
        r_abs_avg_3D[i] = 0;
    }

    // Loop over multiple random walks
    for (int run = 0; run < NUM_RUNS; run++) {
        x_3D[0] = 0;
        y_3D[0] = 0;
        z_3D[0] = 0;
        // Perform a single 3D random walk
        random_walk_3D(NUM_STEPS, DX, DT, D, x_3D, y_3D, z_3D);

        // Accumulate results for each time step
        for (int t = 0; t < NUM_STEPS; t++) {
            double r_squared = pow(x_3D[t], 2) + pow(y_3D[t], 2) + pow(z_3D[t], 2);
            double r_abs = sqrt(r_squared);
            r_squared_avg_3D[t] += r_squared;
            r_abs_avg_3D[t] += r_abs;
        }
    }

    // Compute averages
    for (int t = 0; t < NUM_STEPS; t++) {
        r_squared_avg_3D[t] /= NUM_RUNS;
        r_abs_avg_3D[t] /= NUM_RUNS;
    }

    // Save results to a file
    FILE *file_3D = fopen(filename_3D, "w");


    // Write data to file: time (T), <r�>, <|r|>
    for (int t = 1; t < NUM_STEPS; t++) {
        fprintf(file_3D, "%f %f %f\n", t * dt, r_squared_avg_3D[t], r_abs_avg_3D[t]);
    }

    fclose(file_3D);
    printf("Results saved to %s\n", filename_3D);



    /* This section of the code handles the calculation of averages for self avoiding random walks.
     * two functions write_results_to_file_2D and write_results_to_file_3D are created and called here.
     */


    // Write 2D of self avoiding walks results to file
    write_results_to_file_2D(NUM_STEPS, NUM_RUNS, "sa_2d_results.txt");

    // Write 3D of self avoiding walks results to a file
    write_results_to_file_3D(NUM_STEPS, NUM_STEPS, "sa_3d_results.txt");




    return 0;
}

