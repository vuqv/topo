import sys,os

if len(sys.argv) != 5:
    print('[1] path to cmd file')
    print('[2] path to temp.slurm file')
    print('[3] job name short')
    print('[4] bundel number')
    quit()

if not os.path.exists('./tmp'):
    os.mkdir('./tmp')

temp_slurm = open(sys.argv[2], 'r').readlines()
print(temp_slurm)

tag = sys.argv[3]
bund = int(sys.argv[4])
print(f'bundel size = {bund}')

cmds = open(sys.argv[1], 'r').readlines()
cmds = [c for c in cmds if c != '\n']
print(len(cmds))
cmds = [cmds[i:i+bund] for i in range(0, len(cmds), bund)]
print(cmds)
for cmd_i, cmd in enumerate(cmds):
    print(cmd_i, cmd)
    cmd = ''.join(cmd)
    print(cmd_i, cmd)

    cmd_temp_slurm = temp_slurm.copy()
    cmd_temp_slurm[-1] = f'{cmd}\n'
    cmd_temp_slurm[1] = f'#SBATCH -J c{cmd_i}{tag}\n'

    with open(f'tmp/{cmd_i}{tag}.slurm', 'w') as fh:
        for line in cmd_temp_slurm:
            fh.write(line)
    print(f'SAVED: tmp/{cmd_i}{tag}.slurm -> submitting')
    os.popen(f'sbatch tmp/{cmd_i}{tag}.slurm')

