#!/usr/bin/python

file1=open("tempA.dat","r") #Srecko's
#file1=open("srecko-3mu.txt","r") #Srecko's
list1= file1.readlines()
file1.close()

file2=open("tempB.dat","r") #Cory's
#file2=open("Cory-0e3mu-2012-12-03.dat","r") #Cory's
list2= file2.readlines()
file2.close()

NTT=0
NTF=0
NFT=0

CTT=0
CTF=0
CFT=0

for theline1 in list1:
  rn1=theline1.split(' ')[0]
  lb1=theline1.split(' ')[1]
  id1=theline1.split(' ')[2]
  weight=float(theline1.split(' ')[3])
  #print rn + ' ' + lb + ' ' + id
  found=0
  for theline2 in list2: 
      rn2=theline2.split(' ')[0]
      lb2=theline2.split(' ')[1]
      id2=theline2.split(' ')[2]
      if rn1==rn2 and lb1==lb2 and id1==id2:
        found=1
        break
  if found==0: 
    print 'Cory missing!  '+rn1+ ' '+lb1+ ' '+id1#+','
    NTF += weight
    CTF += 1
  else:
    NTT += weight
    CTT += 1
    
  #print 'Cory missing!  '+theline1 #rn1+ ' '+lb1+ ' '+id1


for theline2 in list2:
  rn2=theline2.split(' ')[0]
  lb2=theline2.split(' ')[1]
  id2=theline2.split(' ')[2]
  weight=float(theline2.split(' ')[3])
  #print rn + ' ' + lb + ' ' + id
  found=0
  for theline1 in list1: 
      rn1=theline1.split(' ')[0]
      lb1=theline1.split(' ')[1]
      id1=theline1.split(' ')[2]
      if rn2==rn1 and lb2==lb1 and id2==id1:
        found=weight
        break
  if found==0:
    print 'Srecko missing!: '+rn2+ ':'+lb2+ ':'+id2#+'\','
    NFT += weight
    CFT += 1
    #print 'Srecko missing!: '+theline2 #rn2+ ' '+lb2+ ' '+id2

print "NTT: "+str(NTT) 
print "NTF: "+str(NTF)
print "NFT: "+str(NFT)

print "CTT: "+str(CTT) 
print "CTF: "+str(CTF)
print "CFT: "+str(CFT)
