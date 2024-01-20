def get_intrinsics():
   
   import os

   intr = []
   for line in open('intrinsics.txt', 'r'):
       for word in line.strip().split(','):
           intr.append(word.strip())

   intr.sort()

   unique = []
   for i in range(len(intr)):
      if i==0: 
          unique.append(intr[i])
      elif intr[i]!=unique[-1]:
          unique.append(intr[i])
    
   print('\n'.join(unique))

get_intrinsics()
