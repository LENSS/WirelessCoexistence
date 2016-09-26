#!/usr/bin/env python
import re
import sys

trans=re.compile('/transmitted_frame_count')
retry=re.compile('/retry_count')

f=open(sys.argv[1], 'r')
w=open('result.txt', 'a')
ls=[]
flag=0
count=0
delimiter=0

for line in f:
  line = line.strip()
  if flag==1:
    #print line
    if count!=0 and count%2!=0:
      ls.append(line)
    flag=0
    if count!=0 and count%2==0:
      ls.append(' '+line)
      s=''.join(ls[count-2:count])
      print s
      w.write(s+'\n')
      delimiter=delimiter+1
      if delimiter==5:
        delimiter=0
        w.write('\n')
  r=trans.search(line)
  if r!=None:
    flag=1
    count=count+1
    #print line
#
#  r=retry.search(line)
#  if r!=None:
#    flag=1
#    #print line
flag=0
count=0
delimiter=0
#
ls[:]=[]
print ''
w.write('\n\n')
f.close()
f=open(sys.argv[1], 'r')
for line in f:
  line = line.strip()
  if flag==1:
    #print line
    if count!=0 and count%2!=0:
      ls.append(line)
    flag=0
    if count!=0 and count%2==0:
      ls.append(' '+line)
      s=''.join(ls[count-2:count])
      print s
      w.write(s+'\n')
      delimiter=delimiter+1
      if delimiter==5:
        delimiter=0
        w.write('\n')
  r=retry.search(line)
  if r!=None:
    flag=1
    count=count+1


w.write('\n================\n\n')
w.close()
f.close()


