import csv,sys
r=list(csv.reader(open(sys.argv[1],encoding='utf-8')))
h=r[0]; bad=[(i,len(x)) for i,x in enumerate(r) if len(x)!=len(h)]
ids=[x[0] for x in r[1:]]
print(len(r)-1,'rows; bad:',bad,'; dup ids:',len(ids)-len(set(ids)))
for x in r[1:]:
    for c in (2,3,4):
        if x[c]: float(x[c])
