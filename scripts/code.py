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
