mpirun -np ${SLURM_NPROCS} pwqmmm.x qmmm.inp 2 |tee out.out
module load qe
#mpirun -np ${SLURM_NPROCS} pw.x <KCl.in|tee out.out
cd qe
mpirun -np ${SLURM_NPROCS} dos.x <dos.in
mpirun -np 30 projwfc.x <prof.in |tee out
ls *Ce*\(f* > f_list
ls *O*\(p* > p_list
sumpdos.x -f f_list |awk 'NR>1{print $1,$2,-$3,$4}'> Ce4f.pdos
sumpdos.x -f p_list |awk 'NR>1{print $1,$2,-$3,$4}'> O2p.pdos
mpirun -np ${SLURM_NPROCS} pp.x < origin.in
mpirun -np ${SLURM_NPROCS} pp.x < pp.in
bader -ref origin.cube chg.cube
~                 
