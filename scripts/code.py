import random

def DNA_RNA_Cod(seq):
  '''
  הפונקציה דואגת שהאותיות תהיינה אחידות (אותיות גדולות) והופכת את רצף ה- DNA המקודד לרצף RNA.
  מקבלת: seq.
  מחזירה: RNA_seq.
  '''
  RNA_seq=""


  if seq[0] != ">":
    seq = seq.upper()
    RNA_seq = seq.replace("T","U")
 
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
