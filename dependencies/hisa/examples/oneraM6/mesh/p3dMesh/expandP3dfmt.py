#!/usr/bin/env python

import re
from shutil import copyfile

    
def expandStr(no, file):
    data  = ""
    lines = file.readlines()
    for line in lines:
        s      = line.strip()
                
        find   = bool(re.search(re.escape(str(no))+'\*', s))

        if (find):
            result = re.findall (re.escape(str(no))+'\*(.*?),', s, re.DOTALL)
            sAdd   = " "
            for i in range(no):
                sAdd += result[0] + "  "
                
            s      = re.sub(re.escape(str(no))+'\*(.*?),', sAdd, s, flags=re.DOTALL)
        
        data   += s+"\n"
    return data


iFile = 'm6wing.x.fmt'
oFile = 'm6wing.p3dfmt'
fact  = [3578,3577,1226,1225,74,25,24,20,8,3,2]

copyfile(iFile, oFile)

for i in fact:
    file  = open(oFile, 'r')
    data  = expandStr(int(i),file)
    text_file = open(oFile, "w")
    text_file.write("%s" % data)
    text_file.close()
    

file  = open(oFile, 'r')
s     = file.read();
s     = re.sub(',', ' ', s, flags=re.DOTALL)
text_file = open(oFile, "w")
text_file.write("%s" % s)
text_file.close()
