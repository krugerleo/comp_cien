

groups = ['CACHE', 'FLOPS_DP', 'MEM']

tamanhos = [1024, 4096, 16384, 65536]


import subprocess

#def itera(*args):
	

for g in groups:
	for t in tamanhos:
		
		cmd = ['likwid-perfctr -m -f -C 0 -g', g, './cgSolver', t, '27 -i 10']
		p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		
		out, err = p.communicate()

