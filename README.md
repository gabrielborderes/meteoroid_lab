#use for monitoring paralel simulation - nvtop

python particles.py



python src/particles_gen_3200Phaenton.py
           
python src/simulate_3200Phaethon.py or mpirun -np 30 python src/simulate_3200Phaethon.py







# Whipple
```bash
python particles.py
mpirun -np 30 python simulate.py
python collect.py
mpirun -np 30 python plot_simulate.py
./compile_plot_simulate.sh data/plots_movie data/plots_movie.mp4
mpirun -np 30 python plot_simulate_close.py
./compile_plot_simulate.sh data/plots_movie_close data/plots_movie_close.mp4
```

# Outburst
```bash
python particles_outburst.py
mpirun -np 30 python simulate_outburst.py
python collect_outburst.py
mpirun -np 30 python plot_simulate_outburst.py
./compile_plot_simulate.sh data/plots_movie_outburst data/plots_movie_outburst.mp4
mpirun -np 30 python plot_simulate_close_outburst.py
./compile_plot_simulate.sh data/plots_movie_close_outburst data/plots_movie_close_outburst.mp4



mpirun -n 30 --map-by :OVERSUBSCRIBE python plot_simulate_test.py 

```
