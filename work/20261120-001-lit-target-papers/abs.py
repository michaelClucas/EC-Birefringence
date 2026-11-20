import re,sys
t=open(sys.argv[1],encoding='utf-8').read()
o=open(sys.argv[2],'w',encoding='utf-8')
for e in t.split('<entry>')[1:]:
    i=re.search(r'<id>(.*?)</id>',e).group(1); ti=re.search(r'<title>(.*?)</title>',e,re.S).group(1); s=re.search(r'<summary>(.*?)</summary>',e,re.S).group(1)
    o.write('## '+i+' '+' '.join(ti.split())+'\n'+' '.join(s.split())+'\n\n')
