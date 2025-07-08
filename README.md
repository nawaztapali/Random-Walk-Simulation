# Random and Self-Avoiding Random Walk Simulation

This repository contains C code written for **Assessment 3 of PHY2027** by *Mohammed Nawaz Tapali* (submitted on 12/12/2024). The program performs simulations of **1D**, **2D**, and **3D random walks**, including **self-avoiding random walks (SAWs)**, on a lattice.

---

## 📌 Project Description

This project simulates the motion of particles on a discrete lattice using:
- **Random Walks (RW)**: A stochastic process where each movement direction is chosen at random.
- **Self-Avoiding Random Walks (SAW)**: A random walk where the particle cannot revisit previously occupied sites.

These are modeled in **1D**, **2D**, and **3D**, allowing analysis of:
- Mean squared displacement ⟨x²⟩
- Mean absolute displacement ⟨|x|⟩
- Visualization of trajectories

---

## 🔧 Features

- Random walks in 1D, 2D, and 3D
- Self-avoiding random walks in 2D and 3D
- Visualization-ready output in text format
- Calculates averages for statistical analysis

---

## 🗂️ Output Files

The code generates the following output files:
- `random_walk_1D.txt`  
- `random_walk_2D.txt`  
- `random_walk_3D.txt`  
- `self_avoiding_random_walk_2D.txt`  
- `self_avoiding_random_walk_3D.txt`  
- `1D_results.txt`, `results_2D.txt`, `results_3D.txt`  
- `sa_2d_results.txt`, `sa_3d_results.txt`

Each file contains position or displacement data to enable analysis or plotting using external tools.

---

## 📁 File Structure

```bash
.
├── random_walk.c            # Main simulation code
├── README.md                # Project description
├── *.txt                    # Output data files
```

---

## 🛠️ Compilation

You can compile the program using `gcc`:

```bash
gcc -o random_walk_sim random_walk.c -lm
./random_walk_sim
```

Make sure you have a C compiler and the math library (`-lm`) available on your system.

---

## 📊 Visualization (Optional)

You can use tools like **Python (matplotlib)**, **Gnuplot**, or **Excel** to visualize the output data from the `.txt` files.

---

## 👤 Author

**Mohammed Nawaz Tapali**  
MSc Physics Student, University of Exeter  
Module: PHY2027 (Computational Physics)

---

## 📜 License

This project is provided for educational and academic use.  
**Do not redistribute, reuse, or submit this work elsewhere without explicit permission from the author.**

All rights reserved © 2024 Mohammed Nawaz Tapali.

---
