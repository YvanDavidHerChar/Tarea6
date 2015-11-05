E := 1
Alpha := 30

DaddysHome.pdf: trayectoria_$(E)_$(Alpha).dat
	python graficando.py trayectoria_$(E)_$(Alpha).dat

trayectoria_$(E)_$(Alpha).dat: particle_in_field.x
	./particle_in_field.x $(E) $(Alpha)

particle_in_field.x: particle.c
	gcc particle.c -lm -o particle_in_field.x


