# --- Import Modules ---
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBList, PDBParser

# Try to import matplotlib (optional)
try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    print("⚠ matplotlib not installed. Skipping visualization step.")

# Setup
Entrez.email = "uttamranjan123@gmail.com"   # Replace with valid email
nuc_accession = "NM_001127511"              # Updated accession number
pdb_id = "1AY7"                             # Updated PDB ID

# --- 1. Parse & Fetch Nucleotide Sequence ---
print("\n--- Fetching Nucleotide Sequence ---")
handle = Entrez.efetch(db="nucleotide", id=nuc_accession, rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()

print("Accession:", record.id)
print("Description:", record.description)

# --- 2. Analyze Features (GC content, CDS extraction) ---
seq = record.seq
gc_content = round((seq.count("G")+seq.count("C"))/len(seq)*100, 2)
print("Total length:", len(seq))
print("GC content:", gc_content, "%")

cds_seq, protein = None, None
for feature in record.features:
    if feature.type == "CDS":
        cds_seq = feature.extract(seq)
        protein = cds_seq.translate(to_stop=True)
        break

if cds_seq:
    print("\n--- Central Dogma ---")
    print("CDS length:", len(cds_seq))
    print("Protein length:", len(protein))
    print("Protein (first 60 aa):", protein[:60])
else:
    print("No CDS found.")
    exit()

# --- 3. Alignment Demo using PairwiseAligner ---
aligner = PairwiseAligner()
aligner.mode = "global"
alns = aligner.align(protein, protein)
print("\n--- Alignment Demo ---")
print("Alignment score:", alns[0].score)
print(alns[0])

# --- 4. Visualization (Optional) ---
if MATPLOTLIB_AVAILABLE:
    bases = ["A","T","G","C"]
    counts = [seq.count(b) for b in bases]
    plt.bar(bases, counts, color="skyblue")
    plt.title("Nucleotide Composition")
    plt.xlabel("Base")
    plt.ylabel("Count")
    plt.savefig("nucleotide_composition.png")
    plt.close()
    print("Plot saved as nucleotide_composition.png")
else:
    print("Visualization skipped (matplotlib not available).")

# --- 5. BLAST (Access biological databases) ---
print("\n--- Running BLASTn (CDS vs nt) ---")
blastn_handle = NCBIWWW.qblast("blastn", "nt", cds_seq)
blastn_record = NCBIXML.read(blastn_handle)
print("Top BLASTn hit:", blastn_record.alignments[0].hit_def)

print("\n--- Running BLASTp (Protein vs nr) ---")
blastp_handle = NCBIWWW.qblast("blastp", "nr", protein)
blastp_record = NCBIXML.read(blastp_handle)
print("Top BLASTp hit:", blastp_record.alignments[0].hit_def)

# --- 6. Fetch PDB Structure ---
print("\n--- Fetching PDB Structure ---")
pdbl = PDBList()
pdb_file = pdbl.retrieve_pdb_file(pdb_id, pdir=".", file_format="pdb")

parser = PDBParser(QUIET=True)
structure = parser.get_structure(pdb_id, pdb_file)
print("PDB structure loaded:", structure)
print("Model count:", len(structure))
print("Chains:", [chain.id for chain in structure[0]])
print("Residues in first chain:", len([res for res in structure[0][list(structure[0].child_dict.keys())[0]]]))

# --- 7. Simple Annotation ---
print("\n--- Annotations (first 5 features) ---")
for i, feature in enumerate(record.features[:5]):
    print(i, feature.type, feature.location)
