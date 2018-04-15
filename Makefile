all: log.system

log.system: system.in
	lmp_mpi -in system.in
