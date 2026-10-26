import random, sys
random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 1)
RATE=float(sys.argv[3]) if len(sys.argv)>3 else 0.002
G=int(sys.argv[2]) if len(sys.argv)>2 else 300000
base=''.join(random.choice('ACGT') for _ in range(G))
def mut(s, rate=RATE):
    out=[];i=0
    while i<len(s):
        r=random.random()
        if r<rate: out.append(random.choice('ACGT')); i+=1
        elif r<rate*1.2: i+=random.randint(1,20)
        elif r<rate*1.4: out.append(''.join(random.choice('ACGT') for _ in range(random.randint(1,20))))
        elif r<rate*1.401: out.append('N'*random.randint(1,300))
        else: out.append(s[i]); i+=1
    return ''.join(out)
with open('collection.txt','w') as f:
    for k in range(12): f.write(mut(base)+'\n')
with open('queries.txt','w') as f:
    for k in range(2): f.write(mut(base,RATE*2)+'\n')
