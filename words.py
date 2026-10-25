#!/usr/bin/env python3
# Turn text files into a sequence of word ids for `hybrid` (32-bit little-endian
# integers, the format `hybrid build` and `hybrid bench` read with width 4).
#
#   python3 words.py OUT.ids [--lower] [--sep] file1.txt [file2.txt ...]
#
# Tokens are maximal runs of letters/digits/apostrophes, and single punctuation
# characters; whitespace is dropped.  Ids are 1, 2, ... in order of decreasing
# frequency.  With --sep, a 0 is written between files.  Also writes OUT.vocab
# (one token per line, line i = id i).
import re, sys, collections, array

args = sys.argv[1:]
if not args: sys.exit(__doc__)
out = args.pop(0)
lower = "--lower" in args; sep = "--sep" in args
files = [a for a in args if a not in ("--lower", "--sep")]
tok = re.compile(r"[\w']+|[^\w\s]", re.UNICODE)
docs = []
for fn in files:
    t = open(fn, encoding="utf-8", errors="replace").read()
    if lower: t = t.lower()
    docs.append(tok.findall(t))
cnt = collections.Counter(w for d in docs for w in d)
vocab = [w for w, _ in cnt.most_common()]
ids = {w: i + 1 for i, w in enumerate(vocab)}
a = array.array("I")
for k, d in enumerate(docs):
    if sep and k: a.append(0)
    a.extend(ids[w] for w in d)
with open(out, "wb") as f: a.tofile(f)
with open(out + ".vocab", "w", encoding="utf-8") as f:
    for w in vocab: f.write(w.replace("\n", " ") + "\n")
print("%d files, %d tokens, %d distinct -> %s" % (len(files), len(a), len(vocab), out))
