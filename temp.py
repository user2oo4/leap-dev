f = open(f'instances_2/sus3.txt', 'w', encoding='utf-8')
f.write('sus3\n')
f.write('10\n')

import random
random.seed(9835)
for i in range(10):
    f.write('0\n')
u = 0
for i in range(10):
    for j in range(i+1,10):
        v = random.randint(5,30)
        if (i>=4 or j<4):
            f.write(f'{i} {j} {-v}\n')
        else:
            f.write(f'{i} {j} {v}\n')
        u+=v
f.write('-1\n')
f.write(f'{-u}\n')

