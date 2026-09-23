"""Fetch SRA run info for PRJNA997998: run accessions, spots, bases, size estimate."""
import re, time, urllib.request

ids = ['28574650','28574649','28574648','28574647','28574646','28574645','28574644','28574643','28574642','28574641','28574640','28574639']
total_bases = 0
for i in ids:
    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={i}&rettype=xml&retmode=xml'
    s = None
    for att in range(4):
        try:
            s = urllib.request.urlopen(url, timeout=60).read().decode('utf-8', errors='ignore')
            break
        except Exception as e:
            time.sleep(3 * (att + 1))
    if s is None:
        print(i, 'FETCH FAILED'); continue
    acc = re.findall(r'accession="(SRR\d+)"', s)
    title = re.findall(r'<TITLE>(.*?)</TITLE>', s)
    spots = re.findall(r'total_spots="(\d+)"', s)
    bases = re.findall(r'total_bases="(\d+)"', s)
    size_est = int(bases[0]) / 2 if bases else 0   # bases/2 ≈ bytes for ~1 byte/base in gz FASTQ
    total_bases += int(bases[0]) if bases else 0
    print(i, acc[0] if acc else '?', title[1][:50] if len(title) > 1 else '?', 'spots:', spots[0] if spots else '?', 'bases(G):', round(int(bases[0])/1e9, 2) if bases else '?', '≈GB(fastq):', round(size_est/1e9, 2))
print('\nTOTAL bases:', round(total_bases / 1e9, 1), 'G | est fastq GB:', round(total_bases / 2 / 1e9, 1))
