$f0=$ARGV[0];
open(F,$f0);
while(<F>){
    ($chr)=split;
    open(SL,">chr${chr}.slum");
    print SL  "#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=4\n#SBATCH --mem=100G\n#SBATCH --time=6:00:00\n#SBATCH --mail-user=xxx\@xxx \n#SBATCH --mail-type=ALL\n#SBATCH -o chr${chr}.log\n";
    print SL "                                \n";
    print SL "ml GCC/8.2.0  CUDA/10.1.105  OpenMPI/3.1.4 R/3.6.0\n";
    print SL "Rscript Trans_mQTL.GTEx.R $chr \n"; # Sort
    close(SL);
    print "${chr}:";
    sleep 1;
    system "sbatch chr${chr}.slum";
}


