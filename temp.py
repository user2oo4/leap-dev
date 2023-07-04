f = open(f'sus.txt', 'w', encoding='utf-8')
f.write('sus\n')
f.write('20\n')

import random
random.seed(123134)
for i in range(20):
    f.write('0\n')
u = 0
for i in range(20):
    for j in range(i+1,20):
        v = random.randint(5,30)
        if (i>=10 or j<10):
            f.write(f'{i} {j} {-v}\n')
        else:
            f.write(f'{i} {j} {v}\n')
        u+=v
f.write('-1\n')
f.write(f'{-u}\n')

