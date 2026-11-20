"""Query the INSPIRE-HEP API, sorted by citations. Usage: python inspire_search.py "<INSPIRE query>" [...]

Network only, no compute. Run via tools/pin.py with PYTHONIOENCODING=utf-8."""
import json,sys,urllib.request,urllib.parse
for q in sys.argv[1:]:
  u="https://inspirehep.net/api/literature?"+urllib.parse.urlencode({'q':q,'size':25,'sort':'mostcited','fields':'titles,authors.full_name,arxiv_eprints,citation_count,publication_info,dois,preprint_date,imprints'})
  d=json.load(urllib.request.urlopen(u))
  print('=== ',q, d['hits']['total'])
  for h in d['hits']['hits']:
    m=h['metadata']; a=m.get('authors',[])
    yr=(m.get('preprint_date') or (m.get('publication_info') or [{}])[0].get('year') or '')
    print(m.get('arxiv_eprints',[{}])[0].get('value','-'), (m.get('dois') or [{}])[0].get('value',''), str(yr)[:4], m.get('citation_count'), '|', m['titles'][0]['title'][:100], '|', ', '.join(x['full_name'] for x in a[:3]))
