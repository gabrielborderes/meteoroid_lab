#use for monitoring paralel simulation - nvtop


# Standard sequence for running the code
```bash

## First step : Prepare particles for simulation 
python src/particles_gen_3200Phaethon.py
           
## Second step: run the simulation with or without mpi           
python src/simulate_3200Phaethon.py 
   or 
mpirun -np 5 python src/simulate_3200Phaethon.py

## Third step: if you runned with mpi   
python src/collect_3200Phaethon.py

## Fourth step: plot the results   
python src/plot_3panels_3200Phaethon.py


```

```
