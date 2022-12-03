all: adg

hw4: adg.cu
	nvcc -o adg adg.cu

clean:
	rm -rf adg
	rm slurm*