"""Verify arXiv IDs and print title/authors/journal/abstract. Usage: python arxiv_verify_ids.py <abstract_chars> <id> [<id> ...]

Network only, no compute. Run via tools/pin.py with PYTHONIOENCODING=utf-8."""
import sys,urllib.request,urllib.parse,time
import xml.etree.ElementTree as ET
ns={'a':'http://www.w3.org/2005/Atom','ax':'http://arxiv.org/schemas/atom'}
ids=sys.argv[2:]; L=int(sys.argv[1])
u="http://export.arxiv.org/api/query?"+urllib.parse.urlencode({'id_list':','.join(ids),'max_results':len(ids)})
r=ET.fromstring(urllib.request.urlopen(u).read())
for e in r.findall('a:entry',ns):
  i=e.find('a:id',ns).text.split('/abs/')[-1]
  t=e.find('a:title',ns)
  if t is None: print('!! NOT FOUND', i); continue
  au=[x.find('a:name',ns).text for x in e.findall('a:author',ns)]
  j=e.find('ax:journal_ref',ns); d=e.find('ax:doi',ns)
  print('##',i, e.find('a:published',ns).text[:4],'|',' '.join(t.text.split()),'|',', '.join(au[:5]), '|', j.text if j is not None else '', d.text if d is not None else '')
  print(' '.join(e.find('a:summary',ns).text.split())[:L]); print()
