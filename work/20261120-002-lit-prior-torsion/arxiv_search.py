"""Query the arXiv API. Usage: python arxiv_search.py "<arXiv search_query>" [...]

Network only, no compute. Run via tools/pin.py with PYTHONIOENCODING=utf-8."""
import sys,urllib.request,urllib.parse,re,time
import xml.etree.ElementTree as ET
ns={'a':'http://www.w3.org/2005/Atom'}
for q in sys.argv[1:]:
  u="http://export.arxiv.org/api/query?"+urllib.parse.urlencode({'search_query':q,'max_results':40,'sortBy':'relevance'})
  r=ET.fromstring(urllib.request.urlopen(u).read())
  print('===',q)
  for e in r.findall('a:entry',ns):
    i=e.find('a:id',ns).text.split('/abs/')[-1]
    t=' '.join(e.find('a:title',ns).text.split())
    au=[x.find('a:name',ns).text for x in e.findall('a:author',ns)]
    print(i, e.find('a:published',ns).text[:4],'|',t[:110],'|',', '.join(au[:3]))
  time.sleep(3)
