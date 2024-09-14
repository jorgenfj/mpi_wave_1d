
# 2.2 Theory Questions

1. **What speed-up/expected speed-up did you get by parallelising the code?**

   The speed-up obtained by parallelising the code depends on several factors, including the problem size, the number of processes used, and the efficiency of communication between processes. In an ideal scenario with no communication overhead and perfect load balancing, the expected speed-up is proportional to the number of processes (linear speed-up). However, due to communication costs and potential load imbalance, the actual speed-up is often less than linear.

   **Sequential:** 5.222492, 4.635370, 4.496746, 5.277175, 4.642154

   **Avg:** 4.8547874

   **MPI_1:** 4.569592, 5.089011, 4.552138, 4.895048, 4.651812

   **Avg:** 4.7515202

   **MPI_4:** 1.675448, 1.702544, 1.666990, 1.675969, 1.788274

   **Avg:** 1.701845

   **MPI_13:** 1.108017, 2.143926, 2.101166, 1.160981, 1.159765

   **Avg:** 1.534771

  Something to note from the provided times:

   The parallel MPI-program run with just 1 process is actually measured as a faster program than the pure sequential one. A reason for this might be because of just using one process there is little overhead incurred with MPI calls, and the sequential program might have been affected by other processes or the operating system itself, causing the additional observed time. 
   
   
   Speedup between MPI with 4 processes and the sequential program:
   $$
   \frac{4.85}{1.70} ≃ 2.85
   $$

   Notice there is a significant speedup, but it is not the ideal linear speedup.

   As the problem is divided into 4 parts, each process performs fewer computations, leading to faster completion. However, the overall performance is still affected by communication overhead and synchronization costs when exchanging boundary data (ghost cells) between processes.

   Speedup between MPI with 13 processes and the sequential program:
   $$
   \frac{4.85}{1.53} ≃ 3.16
   $$

   Here you see that when running with too many (13) processes we are far away from the ideal linear speedup.

   As the number of processes increases, the domain decomposition results in smaller subdomains for each process. This leads to a higher communication-to-computation ratio because the processes exchange boundary data more frequently, and the computation within each subdomain becomes relatively smaller.

   The marginal speedup gained between running with 4 and 13 processes might not be worth it considering the additional resources required running with too many processes.


2. **What is the name for this type of parallelisation (splitting the domain between processes and using ghost cells to communicate)?**

   This type of parallelisation is called **domain decomposition**. In domain decomposition, the computational domain is split into subdomains, each assigned to a different process. Processes use ghost cells (or halo cells) to store boundary data from neighboring subdomains, facilitating communication between processes.

3. **What is the difference between point-to-point communication and collective communication?**

   - **Point-to-point communication** involves data exchange between two specific processes. Examples include `MPI_Send` and `MPI_Recv`. This type of communication is often used when only certain pairs of processes need to exchange data.
   
   - **Collective communication** involves a group of processes within a communicator. Examples include `MPI_Bcast`, `MPI_Gather`, and `MPI_Scatter`. This type of communication is used when data needs to be shared, gathered, or distributed among all processes in a communicator.

4. **Given the following code:**

   ```c
   int* a, b;
   ```
   Which type is `b`?

   `b` is of type `int`. In the declaration `int* a, b;`, only `a` is declared as a pointer to an integer, while `b` is a regular integer.
