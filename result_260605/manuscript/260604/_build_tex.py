#!/usr/bin/env python3
import zipfile, re, os, subprocess, glob, sys

PANDOC = "/usr/local/bin/miniconda3/envs/r_env/bin/pandoc"
STAMP = sys.argv[1]

# 1) extract images image1..7 -> figures_tex/figureN.png
os.makedirs("figures_tex", exist_ok=True)
z = zipfile.ZipFile("Figures.bak.docx")
for n in z.namelist():
    m = re.match(r'word/media/image(\d+)\.(\w+)', n)
    if m:
        idx, ext = int(m.group(1)), m.group(2)
        with open(f"figures_tex/figure{idx}.{ext}", "wb") as fh:
            fh.write(z.read(n))
imgmap = {}
for f in glob.glob("figures_tex/figure*.*"):
    i = int(re.search(r'figure(\d+)', f).group(1))
    imgmap[i] = os.path.basename(f)
print("images:", sorted(imgmap.items()))

# 2) text-body md (sections minus figure legends; strip > notes)
order = ["../sections/00_title_authors.md", "../sections/01_introduction.md",
         "method.md", "result.md", "discussion.md",
         "../sections/05_back_matter.md", "../sections/06_references.md"]
body = []
for p in order:
    for ln in open(p).read().splitlines():
        if ln.strip().startswith('>'):
            continue
        body.append(ln)
    body.append("\n")
open("_body_tex.md", "w").write("\n".join(body))

# 3) pandoc body -> standalone latex
texfile = f"manuscript_{STAMP}.tex"
subprocess.run([PANDOC, "_body_tex.md", "-s", "-t", "latex",
                "-V", "geometry:margin=1in", "-V", "fontsize=11pt",
                "-o", texfile], check=True)
tex = open(texfile).read()

# 4) parse captions -> 7 (heading, body); convert each to latex
caps, head, cb = [], None, []
for ln in open("figure_captions.md").read().splitlines():
    if ln.startswith('## '):
        if head:
            caps.append((head, '\n'.join(cb).strip()))
        head, cb = ln[3:].strip(), []
    elif head is not None and not ln.startswith('>') and not ln.startswith('# '):
        cb.append(ln)
if head:
    caps.append((head, '\n'.join(cb).strip()))

def md2tex_inline(s):
    r = subprocess.run([PANDOC, "-f", "gfm", "-t", "latex"],
                       input=s, capture_output=True, text=True)
    return r.stdout.strip()

blocks = [r"\clearpage", r"\section*{Figures}", ""]
for i, (h, b) in enumerate(caps, start=1):
    capltx = md2tex_inline(f"**{h}**  {b}")
    img = imgmap.get(i, "")
    blocks.append(r"\begin{figure}[htbp]\centering")
    blocks.append(r"\includegraphics[width=\linewidth,height=0.82\textheight,keepaspectratio]{figures_tex/%s}" % img)
    blocks.append(r"\caption*{%s}" % capltx)
    blocks.append(r"\end{figure}")
    blocks.append(r"\clearpage")
    blocks.append("")
figtex = "\n".join(blocks)

# 5) inject packages + figures
pkgs = r"\usepackage{graphicx}" + "\n" + r"\usepackage[skip=4pt]{caption}" + "\n" + r"\usepackage{float}" + "\n"
tex = tex.replace(r"\begin{document}", pkgs + r"\begin{document}", 1)
tex = tex.replace(r"\end{document}", figtex + "\n" + r"\end{document}", 1)
open(texfile, "w").write(tex)
print("\nwrote", texfile, "(", len(tex.split()), "tokens )")
print("figures:", len(caps), "floats with caption*")
