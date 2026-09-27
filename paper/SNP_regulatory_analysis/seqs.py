import re, textwrap
def read_fa(p):
    txt=open(p, encoding='utf-8', errors='ignore').read()
    if txt.lstrip().startswith('{\\rtf'):
        # crude RTF strip
        txt=re.sub(r'\{\\[^{}]*\}','',txt)
        txt=re.sub(r'\\[a-zA-Z]+-?\d* ?','',txt)
        txt=txt.replace('{','').replace('}','')
    out={};name=None
    for line in txt.splitlines():
        line=line.strip()
        if not line: continue
        if line.startswith('>'):
            name=line[1:].split()[0]; out[name]=''
        elif name:
            s=''.join(c for c in line if c.isalpha())
            out[name]+=s
    return out
