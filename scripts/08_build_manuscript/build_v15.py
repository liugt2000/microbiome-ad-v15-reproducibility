# -*- coding: utf-8 -*-
"""Build Manuscript_Microbiome_v15.docx from v15_content.py."""
import re
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

from v15_content import CONTENT, REFS, TABLE1_CAPTION, TABLE1_ROWS, TABLE1_NOTE, FIGS

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, 'results', 'manuscript', 'Manuscript_Microbiome_v15.docx')
os.makedirs(os.path.dirname(OUT), exist_ok=True)
TNR = 'Times New Roman'

# ---------- taxon italic handling ----------
TAXA_PHRASES = [
    'Lachnospiraceae NK4A136 group', 'Ruminococcus gnavus group', 'Ruminococcus torques group',
    'Ruminococcus gauvreauii group', 'Eubacterium hallii group', 'Eubacterium eligens group',
    'Eubacterium ventriosum group', 'Eubacterium nodatum group', 'Christensenellaceae R-7 group',
    'Family XIII AD3011 group', 'Clostridium sensu stricto 1', 'Ruminococcaceae UCG013',
    'Lachnospiraceae UCG010', 'NK4A214 group', 'Escherichia\u2013Shigella',
    'Faecalibacterium prausnitzii', 'Roseburia intestinalis', 'Eubacterium rectale',
    'Akkermansia muciniphila', 'Clostridium butyricum', 'Incertae Sedis',
    'Faecalibacterium', 'Roseburia', 'Blautia', 'Butyricicoccus', 'Coprococcus',
    'Agathobacter', 'Subdoligranulum', 'Lachnospira', 'Dorea', 'Fusicatenibacter',
    'Streptococcus', 'Veillonella', 'Rothia', 'Actinomyces', 'Granulicatella',
    'Haemophilus', 'Marvinbryantia', 'Monoglobus', 'Anaerostipes', 'Bacteroides',
    'Alistipes', 'Akkermansia', 'Bilophila', 'Desulfovibrio', 'Prevotella',
    'Collinsella', 'Fusobacterium', 'Intestinimonas', 'Lactobacillus',
    'Bifidobacterium', 'Enterococcus', 'Dialister', 'Romboutsia', 'Lachnoclostridium',
    'Flavonifractor', 'Holdemanella', 'Senegalimassilia', 'Methanobrevibacter',
    'Oscillibacter', 'Clostridium', 'Enterobacter', 'Klebsiella', 'Hungatella',
    'Eggerthella', 'Anaerotruncus', 'Odoribacter', 'Pseudobutyrivibrio',
    'Butyrivibrio', 'Christensenella', 'Eubacterium', 'Ruminococcus',
    'Lachnospiraceae', 'Ruminococcaceae', 'Christensenellaceae', 'Oscillospiraceae',
    'Oscillospirales', 'Enterobacteriaceae', 'Acetobacterium', 'UCG-005',
]
TAXA_PHRASES.sort(key=len, reverse=True)
_TAXA_RE = re.compile('(' + '|'.join(re.escape(t) for t in TAXA_PHRASES) + ')')

def add_rich(p, text, size=None, bold=False):
    """Add text to paragraph p, italicizing taxon names."""
    parts = _TAXA_RE.split(text)
    for part in parts:
        if not part:
            continue
        r = p.add_run(part)
        r.italic = bool(_TAXA_RE.fullmatch(part))
        if size:
            r.font.size = Pt(size)
        if bold:
            r.font.bold = True

# ---------- document setup ----------
doc = Document()

sec = doc.sections[0]
sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
sec.top_margin = sec.bottom_margin = sec.left_margin = sec.right_margin = Cm(2.5)

def _font(style, name=TNR, size=12, bold=None, color=(0, 0, 0)):
    style.font.name = name
    style.font.size = Pt(size)
    if bold is not None:
        style.font.bold = bold
    style.font.color.rgb = RGBColor(*color)
    rpr = style.element.get_or_add_rPr()
    rf = rpr.find(qn('w:rFonts'))
    if rf is None:
        rf = OxmlElement('w:rFonts')
        rpr.append(rf)
    rf.set(qn('w:ascii'), name)
    rf.set(qn('w:hAnsi'), name)
    rf.set(qn('w:eastAsia'), name)

st = doc.styles['Normal']
_font(st, size=12)
st.paragraph_format.line_spacing = 1.5
st.paragraph_format.space_after = Pt(8)
st.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY

h1 = doc.styles['Heading 1']
_font(h1, size=14, bold=True)
h1.paragraph_format.space_before = Pt(18)
h1.paragraph_format.space_after = Pt(8)
h1.paragraph_format.keep_with_next = True

h2 = doc.styles['Heading 2']
_font(h2, size=12, bold=True)
h2.paragraph_format.space_before = Pt(10)
h2.paragraph_format.space_after = Pt(4)
h2.paragraph_format.keep_with_next = True

# footer page number
footer_p = sec.footer.paragraphs[0]
footer_p.alignment = WD_ALIGN_PARAGRAPH.CENTER
run = footer_p.add_run()
begin = OxmlElement('w:fldChar')
begin.set(qn('w:fldCharType'), 'begin')
instr = OxmlElement('w:instrText')
instr.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
instr.text = ' PAGE '
end = OxmlElement('w:fldChar')
end.set(qn('w:fldCharType'), 'end')
for e in (begin, instr, end):
    run._r.append(e)
run.font.name = TNR
run.font.size = Pt(10)

# title
t = doc.add_paragraph()
t.alignment = WD_ALIGN_PARAGRAPH.CENTER
t.paragraph_format.space_after = Pt(18)
t.paragraph_format.line_spacing = 1.3
tr = t.add_run(CONTENT[0][1])
tr.font.name = TNR
tr.font.size = Pt(16)
tr.font.bold = True

# ---------- body ----------
refs_started = False
for kind, payload in CONTENT[1:]:
    if kind == 'h1':
        refs_started = payload.startswith('11.')
        doc.add_paragraph(payload, style='Heading 1')
    elif kind == 'h2':
        doc.add_paragraph(payload, style='Heading 2')
    elif kind == 'p':
        p = doc.add_paragraph()
        add_rich(p, payload)
    elif kind == 'pb':
        label, rest = payload
        p = doc.add_paragraph()
        r = p.add_run(label)
        r.font.bold = True
        add_rich(p, rest)

# ---------- References ----------
for i, ref in enumerate(REFS, 1):
    p = doc.add_paragraph()
    pf = p.paragraph_format
    pf.left_indent = Cm(0.75)
    pf.first_line_indent = Cm(-0.75)
    pf.line_spacing = 1.15
    pf.space_after = Pt(3)
    add_rich(p, '[%d] %s' % (i, ref), size=10)

# ---------- Table 1 ----------
p = doc.add_paragraph()
p.paragraph_format.keep_with_next = True
p.paragraph_format.space_before = Pt(12)
r = p.add_run('Table 1. ')
r.font.bold = True
r.font.size = Pt(10.5)
add_rich(p, TABLE1_CAPTION, size=10.5)

# ---------- Table 1 ----------
table = doc.add_table(rows=len(TABLE1_ROWS), cols=5)
table.style = 'Table Grid'
table.alignment = WD_TABLE_ALIGNMENT.CENTER
widths = [Cm(3.3), Cm(4.9), Cm(3.8), Cm(2.1), Cm(1.9)]

def shade(cell, hexcolor):
    tcPr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement('w:shd')
    shd.set(qn('w:val'), 'clear')
    shd.set(qn('w:fill'), hexcolor)
    tcPr.append(shd)

for ri, row in enumerate(TABLE1_ROWS):
    tr = table.rows[ri]
    tr.height = Cm(0.6)
    # keep row on one page
    trPr = tr._tr.get_or_add_trPr()
    cant = OxmlElement('w:cantSplit')
    trPr.append(cant)
    if ri == 0:
        th = OxmlElement('w:tblHeader')
        trPr.append(th)
    for ci, val in enumerate(row):
        cell = tr.cells[ci]
        cell.width = widths[ci]
        p = cell.paragraphs[0]
        p.paragraph_format.line_spacing = 1.0
        p.paragraph_format.space_after = Pt(2)
        p.paragraph_format.first_line_indent = Cm(0)
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER if (ri == 0 or ci in (2, 3, 4)) else WD_ALIGN_PARAGRAPH.LEFT
        add_rich(p, val, size=8.5, bold=(ri == 0))
        if ri == 0:
            shade(cell, 'F2F2F2')

note = doc.add_paragraph()
note.paragraph_format.space_before = Pt(6)
note.paragraph_format.line_spacing = 1.15
add_rich(note, TABLE1_NOTE, size=9)

# ---------- Figure legends ----------
doc.add_paragraph('12. Figure Legends', style='Heading 1')
for label, text in FIGS:
    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(10)
    p.paragraph_format.line_spacing = 1.3
    r = p.add_run(label)
    r.font.bold = True
    add_rich(p, text)

# document properties
cp = doc.core_properties
cp.title = 'Integrative Multi-Omics and Multi-Cohort Meta-Analysis Reveal Depleted Butyrate- and Hydrogen-Metabolizing Microbial Consortia and Reduced Fecal Taurine Across the Alzheimer\u2019s Disease Continuum'
cp.author = '[Author list to be completed]'
cp.subject = 'Manuscript v15'

doc.save(OUT)
print('SAVED:', OUT)
print('paragraphs:', len(doc.paragraphs), 'tables:', len(doc.tables))
