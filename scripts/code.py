import random

def DNA_RNA_Cod(seq):
  '''
  הפונקציה דואגת שהאותיות תהיינה אחידות (אותיות גדולות) והופכת את רצף ה- DNA המקודד לרצף RNA.
  מקבלת: seq.
  מחזירה: RNA_seq.
  '''
  RNA_seq=""

  for line in seq:
    line = line.upper()
    line = line.rstrip('\r\n')
    
    if line[0] != ">":
      rna_line = line.replace("T","U")
      RNA_seq += rna_line
    
  return RNA_seq
#------------------------------------------------
 
def Read_dict(fl):
  '''
  הפונקציה קוראת לתוךdictionary  את המיפוי בין הקודונים לחומצות אמינו מהקובץ.
  מקבלת: fl.
  מחזירה: RNA_codon_table.
  '''
  global RNA_codon_table

  for line in fl:
    line = line.rstrip('\r\n')
    line_list = line.split()
    
    codon_key = line_list[0]
    amino_acids_value = line_list[1]
    
    RNA_codon_table[codon_key] = amino_acids_value
    
  return RNA_codon_table
#------------------------------------------------
  
def RNA_prot(seq):
  '''
  הפונקציה מתרגמת את רצף ה- RNA לרצף חלבון.
  מקבלת: seq.
  מחזירה: protein_seq.
  '''
  protein_seq = ""
  start = False
  for i in range(0, len(seq), 3):
    codon = seq[i:i+3]
    if len(codon) < 3:
      break
    
    curr_amino_acid = RNA_codon_table.get(codon,"-")
    
    if curr_amino_acid == "M":
      start = True
      
    if start:
      if curr_amino_acid == "*":
        protein_seq = protein_seq + "*"
        break
      
      protein_seq += curr_amino_acid
  
  return protein_seq
#------------------------------------------------

def Insert_DNA(seq):
  '''
  הפונקציה תכניס במיקום אקראי לרצף ה- DNA של הגן נוקלאוטיד נוסף.
  מקבלת: seq.
  מחזירה: change_genome.
  '''
  nucleotide_list = ['T','G','C','A']
  
  rand_nucleotide = random.choice(nucleotide_list)
  rand_num = random.randrange(0,len(seq))
  
  change_genome = seq[0:rand_num]+ rand_nucleotide + seq[rand_num:]
  
  return change_genome
#------------------------------------------------
  
def Delete_DNA(seq):
  '''
  הפונקציה תחסיר נוקלאוטיד במיקום רנדומאלי.
  מקבלת: seq.
  מחזירה: change_genome.
  '''
  rand_num = random.randrange(0,len(seq))
  rand_nucleotide = seq[rand_num]
  
  change_genome = seq[:rand_num] + seq[rand_num + 1:]

  return change_genome
#------------------------------------------------

def Mutate_DNA(seq):
  '''
  הפונקציה בוחרת מיקום אקראי ברצף גנום ה- HIV ומחליפה במיקום זה את הנוקלאוטיד באופן רנדומלי ל- A/C/G/T.
  מקבלת: seq.
  מחזירה:change_genome
  '''
  nucleotide_list = ['T','G','C','A']
  
  rand_nucleotide = random.choice(nucleotide_list)
  rand_num = random.randrange(0,len(seq))
  
  if seq[rand_num] != rand_nucleotide:
    change_genome = seq[0:rand_num]+ rand_nucleotide + seq[(rand_num+1):]
  
  else:
    nucleotide_list.remove(rand_nucleotide)
    rand_nucleotide = random.choice(nucleotide_list)
    change_genome = seq[0:rand_num]+ rand_nucleotide + seq[(rand_num+1):]
  return change_genome 

# תוכנית ראשית.
global RNA_codon_table
RNA_codon_table = {}

# פתיחת הקבצים
p53_seq = open('data/p53_sequence.fa', 'r')
codon_file = open('data/codon_AA.txt', 'r')

# קריאה לפונקציה
Read_dict(codon_file)

p53_genome = ""
new_HIV_genome = ""
# קריאת הקובץ
for line in p53_seq:
  line = line.rstrip('\r\n')
  if line == "":
    continue
  # רצף ה DNA מופיע בשורות שאינן מתחילות בסימן "<" לכן "נדלג" על שורה זו
  if line[0] == ">":
    continue
  
  p53_genome = p53_genome + line

print(p53_genome[2:20]) 
gene_as_RNA = DNA_RNA_Cod(p53_genome)
protein = RNA_prot(gene_as_RNA)
print(protein)      
print(len(p53_genome))

#
for i in range(3):
  num = random.randrange(3)
  if num == 1:
    p53_genome = Mutate_DNA(p53_genome)
    print(1)
  elif num == 2:
    p53_genome = Insert_DNA(p53_genome)
    print(2)
  else:
    p53_genome = Delete_DNA(p53_genome)
    print(3)


print(len(p53_genome))
