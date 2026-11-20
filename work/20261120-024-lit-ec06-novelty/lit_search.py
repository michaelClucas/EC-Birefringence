"""arXiv / INSPIRE queries for tasks 013 and 024 (copied from 013). Network only, no compute.

Usage:
  python lit_search.py arxiv "<arXiv search_query>" [...]
  python lit_search.py inspire "<INSPIRE query>" [...]     (sorted most recent)
  python lit_search.py ids <arXiv id> [...]                 (verify id -> title, authors, date)
Run via tools/pin.py with PYTHONIOENCODING=utf-8.
"""
import json, sys, time, urllib.parse, urllib.request
import xml.etree.ElementTree as ET

NS = {'a': 'http://www.w3.org/2005/Atom'}
UA = {'User-Agent': 'lit-scout/0.1'}


def get(u):
    return urllib.request.urlopen(urllib.request.Request(u, headers=UA), timeout=60).read()


def arxiv(params, abstracts=False):
    r = ET.fromstring(get("http://export.arxiv.org/api/query?" + urllib.parse.urlencode(params)))
    for e in r.findall('a:entry', NS):
        i = e.find('a:id', NS).text.split('/abs/')[-1]
        t = ' '.join(e.find('a:title', NS).text.split())
        au = [x.find('a:name', NS).text for x in e.findall('a:author', NS)]
        print(i, e.find('a:published', NS).text[:10], '|', t[:120], '|', ', '.join(au[:3]))
        if abstracts:
            print('   ', ' '.join(e.find('a:summary', NS).text.split()))
            print()
    time.sleep(3)


def inspire(q):
    u = "https://inspirehep.net/api/literature?" + urllib.parse.urlencode({
        'q': q, 'size': 60, 'sort': 'mostrecent',
        'fields': 'titles,authors.full_name,arxiv_eprints,citation_count,preprint_date'})
    d = json.loads(get(u))
    print('=== ', q, d['hits']['total'])
    for h in d['hits']['hits']:
        m = h['metadata']; a = m.get('authors', [])
        print(m.get('arxiv_eprints', [{}])[0].get('value', '-'), str(m.get('preprint_date', ''))[:7],
              m.get('citation_count'), '|', m['titles'][0]['title'][:110], '|',
              ', '.join(x['full_name'] for x in a[:3]))
    time.sleep(1)


mode, args = sys.argv[1], sys.argv[2:]
if mode == 'arxiv':
    for q in args:
        print('===', q)
        arxiv({'search_query': q, 'max_results': 30, 'sortBy': 'relevance'})
elif mode == 'inspire':
    for q in args:
        inspire(q)
elif mode == 'ids':
    arxiv({'id_list': ','.join(args), 'max_results': len(args)}, abstracts=True)
elif mode == 'inspire_id':
    # INSPIRE record for an arXiv id, then all papers citing it
    for a in args:
        d = json.loads(get("https://inspirehep.net/api/arxiv/" + a + "?fields=control_number,titles,citation_count"))
        rec = d['metadata']['control_number']
        print('### ', a, 'record', rec, 'citations', d['metadata'].get('citation_count'))
        inspire('refersto:recid:%s' % rec)
elif mode == 'citeq':
    # citeq <arXiv id> "<extra INSPIRE query or ->"... : papers citing the id, optionally filtered
    a, extras = args[0], args[1:] or ['-']
    d = json.loads(get("https://inspirehep.net/api/arxiv/" + a + "?fields=control_number"))
    rec = d['metadata']['control_number']
    for x in extras:
        inspire('refersto:recid:%s' % rec + ('' if x == '-' else ' and ' + x))
